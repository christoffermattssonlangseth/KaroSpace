"""Pathway enrichment helpers for compact viewer exports."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.stats import hypergeom

from .pseudobulk import _adjust_pvalues, _json_float


DEFAULT_REACTOME_LIBRARY = "Reactome_2022"


def _clean_gene(value: Any) -> str:
    return str(value or "").strip()


def _finite_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def _dedupe(values: Iterable[str]) -> List[str]:
    seen = set()
    result: List[str] = []
    for value in values:
        token = _clean_gene(value)
        if not token:
            continue
        key = token.lower()
        if key in seen:
            continue
        seen.add(key)
        result.append(token)
    return result


def _parse_gmt_file(path: Union[str, Path]) -> Dict[str, List[str]]:
    file_path = Path(path).expanduser()
    gene_sets: Dict[str, List[str]] = {}
    with file_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0].strip()
            if not term:
                continue
            genes = _dedupe(parts[2:])
            if genes:
                gene_sets[term] = genes
    return gene_sets


def _load_reactome_from_gseapy(organism: str = "Human") -> Tuple[Dict[str, List[str]], str]:
    try:
        import gseapy as gp  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on optional dependency
        raise RuntimeError(
            "Reactome pathway enrichment requires gseapy when --pathway-gmt is not supplied. "
            "Install dependencies or pass a local Reactome GMT with --pathway-gmt."
        ) from exc

    library_name = DEFAULT_REACTOME_LIBRARY
    try:
        library_names = list(gp.get_library_name(organism=organism))
        reactome_names = [
            str(name) for name in library_names if "reactome" in str(name).lower()
        ]
        if reactome_names:
            library_name = sorted(
                reactome_names,
                key=lambda name: (
                    max([int(part) for part in "".join(ch if ch.isdigit() else " " for ch in name).split()] or [0]),
                    name,
                ),
                reverse=True,
            )[0]
    except Exception:
        library_name = DEFAULT_REACTOME_LIBRARY

    try:
        raw_sets = gp.get_library(name=library_name, organism=organism)
    except TypeError:
        raw_sets = gp.get_library(library_name, organism=organism)
    gene_sets = {
        str(term): _dedupe(genes)
        for term, genes in dict(raw_sets).items()
        if _dedupe(genes)
    }
    if not gene_sets:
        raise RuntimeError(f"No genes were loaded from Reactome library '{library_name}'.")
    return gene_sets, library_name


def load_pathway_gene_sets(
    pathway_gmt: Optional[Union[str, Sequence[str]]] = None,
    *,
    organism: str = "Human",
) -> Tuple[Dict[str, List[str]], Dict[str, Any]]:
    """Load pathway gene sets from GMT files or the default Reactome library."""
    if pathway_gmt is None or pathway_gmt == "":
        gene_sets, source_name = _load_reactome_from_gseapy(organism=organism)
        return gene_sets, {
            "source": "reactome",
            "library": source_name,
            "organism": organism,
            "gene_set_count": len(gene_sets),
        }

    paths = [pathway_gmt] if isinstance(pathway_gmt, str) else list(pathway_gmt)
    combined: Dict[str, List[str]] = {}
    loaded_files: List[str] = []
    for raw_path in paths:
        token = str(raw_path or "").strip()
        if not token:
            continue
        gene_sets = _parse_gmt_file(token)
        prefix = Path(token).expanduser().stem
        for term, genes in gene_sets.items():
            key = term if term not in combined else f"{prefix}:{term}"
            combined[key] = genes
        loaded_files.append(str(Path(token).expanduser()))
    if not combined:
        raise RuntimeError("No pathway gene sets were loaded from --pathway-gmt.")
    return combined, {
        "source": "gmt",
        "files": loaded_files,
        "gene_set_count": len(combined),
    }


def _match_gene_sets(
    gene_sets: Mapping[str, Sequence[str]],
    universe_genes: Sequence[str],
    *,
    min_overlap: int,
    max_size: int,
) -> Dict[str, List[str]]:
    universe = [_clean_gene(gene) for gene in universe_genes if _clean_gene(gene)]
    lower_to_gene = {gene.lower(): gene for gene in universe}
    matched: Dict[str, List[str]] = {}
    for term, genes in gene_sets.items():
        values = []
        seen = set()
        for gene in genes:
            raw = _clean_gene(gene)
            if not raw:
                continue
            match = lower_to_gene.get(raw.lower())
            if not match or match in seen:
                continue
            seen.add(match)
            values.append(match)
        if len(values) >= int(min_overlap) and len(values) <= int(max_size):
            matched[str(term)] = values
    return matched


def _result_gene_table(result: Mapping[str, Any]) -> List[Dict[str, Any]]:
    genes = result.get("genes") if isinstance(result.get("genes"), list) else []
    log2fc = result.get("log2foldchanges")
    if not isinstance(log2fc, list):
        log2fc = result.get("logfoldchanges") if isinstance(result.get("logfoldchanges"), list) else []
    pvals_adj = result.get("pvals_adj") if isinstance(result.get("pvals_adj"), list) else []
    pvals = result.get("pvals") if isinstance(result.get("pvals"), list) else []
    scores = result.get("scores") if isinstance(result.get("scores"), list) else []
    pct_source = result.get("pct_source") if isinstance(result.get("pct_source"), list) else []
    pct_reference = result.get("pct_reference") if isinstance(result.get("pct_reference"), list) else []

    rows: List[Dict[str, Any]] = []
    for idx, gene in enumerate(genes):
        token = _clean_gene(gene)
        if not token:
            continue
        fc = _finite_float(log2fc[idx] if idx < len(log2fc) else None)
        padj = _finite_float(pvals_adj[idx] if idx < len(pvals_adj) else None)
        pval = _finite_float(pvals[idx] if idx < len(pvals) else None)
        score = _finite_float(scores[idx] if idx < len(scores) else None)
        if score is None:
            p_for_rank = padj if padj is not None and padj > 0 else pval
            if fc is not None and p_for_rank is not None and p_for_rank > 0:
                score = float(np.sign(fc) * -np.log10(max(p_for_rank, np.nextafter(0.0, 1.0))))
            elif fc is not None:
                score = fc
        source_pct = _finite_float(pct_source[idx] if idx < len(pct_source) else None)
        reference_pct = _finite_float(pct_reference[idx] if idx < len(pct_reference) else None)
        rows.append({
            "gene": token,
            "log2fc": fc,
            "padj": padj,
            "pval": pval,
            "score": score,
            "pct_source": source_pct,
            "pct_reference": reference_pct,
        })
    return rows


def _ora_rows(
    rows: Sequence[Mapping[str, Any]],
    matched_sets: Mapping[str, Sequence[str]],
    *,
    direction: str,
    padj_cutoff: float,
    log2fc_cutoff: float,
    min_pct_expressed: float,
    min_overlap: int,
    top_n: int,
) -> List[Dict[str, Any]]:
    universe = {str(row["gene"]) for row in rows}
    def _passes_pct(row: Mapping[str, Any]) -> bool:
        if min_pct_expressed <= 0:
            return True
        source_pct = _finite_float(row.get("pct_source")) or 0.0
        reference_pct = _finite_float(row.get("pct_reference")) or 0.0
        return source_pct >= min_pct_expressed or reference_pct >= min_pct_expressed

    if direction == "up":
        query = {
            str(row["gene"])
            for row in rows
            if row.get("padj") is not None
            and row.get("log2fc") is not None
            and _passes_pct(row)
            and float(row["padj"]) < padj_cutoff
            and float(row["log2fc"]) >= log2fc_cutoff
        }
    else:
        query = {
            str(row["gene"])
            for row in rows
            if row.get("padj") is not None
            and row.get("log2fc") is not None
            and _passes_pct(row)
            and float(row["padj"]) < padj_cutoff
            and float(row["log2fc"]) <= -log2fc_cutoff
        }
    if len(universe) < 1 or len(query) < int(min_overlap):
        return []

    pvals: List[float] = []
    raw: List[Dict[str, Any]] = []
    M = len(universe)
    n = len(query)
    for term, genes in matched_sets.items():
        pathway_genes = set(genes)
        K = len(pathway_genes)
        overlap_genes = sorted(query & pathway_genes)
        k = len(overlap_genes)
        if k < int(min_overlap):
            continue
        pval = float(hypergeom.sf(k - 1, M, K, n))
        pvals.append(pval)
        background_non_pathway = max(M - K, 0)
        query_non_overlap = max(n - k, 0)
        odds_ratio = ((k + 0.5) * (background_non_pathway - query_non_overlap + 0.5)) / (
            (query_non_overlap + 0.5) * (K - k + 0.5)
        )
        raw.append({
            "term": str(term),
            "direction": direction,
            "pval": pval,
            "odds_ratio": odds_ratio,
            "overlap": k,
            "query_size": n,
            "pathway_size": K,
            "genes": overlap_genes,
        })

    if not raw:
        return []
    adjusted = _adjust_pvalues(np.asarray(pvals, dtype=float), "fdr_bh")
    for entry, padj in zip(raw, adjusted):
        entry["padj"] = float(padj) if np.isfinite(padj) else None
    raw.sort(key=lambda item: (
        float(item["padj"]) if item.get("padj") is not None else float("inf"),
        float(item["pval"]),
        -int(item["overlap"]),
        str(item["term"]),
    ))
    return [_jsonify_pathway_row(entry) for entry in raw[: int(top_n)]]


def _running_enrichment_score(scores: np.ndarray, hit_mask: np.ndarray) -> Tuple[float, np.ndarray]:
    n_genes = int(scores.size)
    n_hits = int(hit_mask.sum())
    if n_genes == 0 or n_hits == 0 or n_hits == n_genes:
        return 0.0, np.zeros(n_genes, dtype=float)
    weights = np.abs(scores)
    hit_weight_sum = float(weights[hit_mask].sum())
    if hit_weight_sum <= 0:
        weights = np.ones_like(scores)
        hit_weight_sum = float(n_hits)
    running = np.where(hit_mask, weights / hit_weight_sum, -1.0 / float(n_genes - n_hits))
    cumulative = np.cumsum(running)
    max_es = float(cumulative.max())
    min_es = float(cumulative.min())
    return (max_es if abs(max_es) >= abs(min_es) else min_es), cumulative


def _sample_indices(length: int, max_points: int, extra: Iterable[int] = ()) -> List[int]:
    if length <= 0:
        return []
    limit = max(2, int(max_points))
    if length <= limit:
        indices = set(range(length))
    else:
        indices = set(int(round(value)) for value in np.linspace(0, length - 1, limit))
    for idx in extra:
        if 0 <= int(idx) < length:
            indices.add(int(idx))
    return sorted(indices)


def _gsea_plot_payload(
    scores: np.ndarray,
    running: np.ndarray,
    hit_indices: Sequence[int],
    peak_idx: int,
) -> Dict[str, Any]:
    """Return compact data for a classic GSEA enrichment-profile plot."""
    n_genes = int(len(scores))
    profile_indices = _sample_indices(n_genes, 180, extra=[0, peak_idx, n_genes - 1])
    metric_indices = _sample_indices(n_genes, 100, extra=[0, n_genes - 1])
    non_positive = np.flatnonzero(scores <= 0)
    zero_cross = int(non_positive[0]) if non_positive.size else None
    return {
        "rank_count": n_genes,
        "peak_rank": int(peak_idx),
        "zero_cross_rank": zero_cross,
        "hit_indices": [int(idx) for idx in hit_indices],
        "running_profile": [
            [int(idx), _json_float(float(running[idx]), 6)]
            for idx in profile_indices
        ],
        "rank_metric_profile": [
            [int(idx), _json_float(float(scores[idx]), 6)]
            for idx in metric_indices
        ],
    }


def _gsea_rows(
    rows: Sequence[Mapping[str, Any]],
    matched_sets: Mapping[str, Sequence[str]],
    *,
    min_overlap: int,
    top_n: int,
    permutations: int,
    seed: int,
) -> Dict[str, List[Dict[str, Any]]]:
    ranked = [
        (str(row["gene"]), float(row["score"]))
        for row in rows
        if row.get("score") is not None and np.isfinite(float(row["score"]))
    ]
    ranked.sort(key=lambda item: (-item[1], item[0]))
    if len(ranked) < max(2, int(min_overlap)):
        return {"positive": [], "negative": []}

    genes = [gene for gene, _score in ranked]
    gene_index = {gene: idx for idx, gene in enumerate(genes)}
    scores = np.asarray([score for _gene, score in ranked], dtype=float)
    rng = np.random.default_rng(int(seed))
    raw: List[Dict[str, Any]] = []
    pvals: List[float] = []
    n_perm = max(0, int(permutations))

    for term, pathway_genes in matched_sets.items():
        hit_indices = sorted(gene_index[gene] for gene in pathway_genes if gene in gene_index)
        if len(hit_indices) < int(min_overlap) or len(hit_indices) >= len(genes):
            continue
        hit_mask = np.zeros(len(genes), dtype=bool)
        hit_mask[hit_indices] = True
        es, running = _running_enrichment_score(scores, hit_mask)
        if es == 0:
            continue
        if n_perm > 0:
            null_scores = np.empty(n_perm, dtype=float)
            for perm_idx in range(n_perm):
                random_hits = rng.choice(len(genes), size=len(hit_indices), replace=False)
                random_mask = np.zeros(len(genes), dtype=bool)
                random_mask[random_hits] = True
                null_scores[perm_idx], _ = _running_enrichment_score(scores, random_mask)
            if es >= 0:
                denom_values = null_scores[null_scores > 0]
                denom = float(np.mean(denom_values)) if denom_values.size else float(np.mean(np.abs(null_scores)))
                pval = float((np.count_nonzero(null_scores >= es) + 1) / (n_perm + 1))
            else:
                denom_values = np.abs(null_scores[null_scores < 0])
                denom = float(np.mean(denom_values)) if denom_values.size else float(np.mean(np.abs(null_scores)))
                pval = float((np.count_nonzero(null_scores <= es) + 1) / (n_perm + 1))
            nes = float(es / denom) if denom > 0 else float(es)
        else:
            pval = 1.0
            nes = float(es)
        peak_idx = int(np.argmax(running) if es >= 0 else np.argmin(running))
        if es >= 0:
            leading = [gene for idx, gene in enumerate(genes[: peak_idx + 1]) if hit_mask[idx]]
            direction = "positive"
        else:
            leading = [gene for idx, gene in enumerate(genes[peak_idx:]) if hit_mask[peak_idx + idx]]
            direction = "negative"
        raw.append({
            "term": str(term),
            "direction": direction,
            "es": es,
            "nes": nes,
            "pval": pval,
            "overlap": len(hit_indices),
            "pathway_size": len(hit_indices),
            "leading_edge": leading[:50],
            "_hit_indices": hit_indices,
            "_peak_idx": peak_idx,
        })
        pvals.append(pval)

    if not raw:
        return {"positive": [], "negative": []}
    adjusted = _adjust_pvalues(np.asarray(pvals, dtype=float), "fdr_bh")
    for entry, padj in zip(raw, adjusted):
        entry["padj"] = float(padj) if np.isfinite(padj) else None
    raw.sort(key=lambda item: (
        float(item["padj"]) if item.get("padj") is not None else float("inf"),
        -abs(float(item["nes"])),
        str(item["term"]),
    ))
    positive_entries = [entry for entry in raw if entry["direction"] == "positive"][: int(top_n)]
    negative_entries = [entry for entry in raw if entry["direction"] == "negative"][: int(top_n)]
    for entry in positive_entries + negative_entries:
        hit_indices = entry.pop("_hit_indices", [])
        peak_idx = int(entry.pop("_peak_idx", 0))
        hit_mask = np.zeros(len(scores), dtype=bool)
        hit_mask[list(hit_indices)] = True
        _es, running = _running_enrichment_score(scores, hit_mask)
        entry.update(_gsea_plot_payload(scores, running, hit_indices, peak_idx))
    positive = [_jsonify_pathway_row(entry) for entry in positive_entries]
    negative = [_jsonify_pathway_row(entry) for entry in negative_entries]
    return {"positive": positive, "negative": negative}


def _jsonify_pathway_value(value: Any) -> Any:
    if isinstance(value, (list, tuple)):
        return [_jsonify_pathway_value(item) for item in value]
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return _json_float(float(value), 6)
    return str(value) if value is not None else None


def _jsonify_pathway_row(entry: Mapping[str, Any]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in entry.items():
        result[key] = _jsonify_pathway_value(value)
    return result


def enrich_pseudobulk_result(
    result: Mapping[str, Any],
    gene_sets: Mapping[str, Sequence[str]],
    *,
    top_n: int = 20,
    min_overlap: int = 3,
    max_pathway_size: int = 500,
    gsea_permutations: int = 100,
    seed: int = 0,
) -> Optional[Dict[str, Any]]:
    """Compute compact ORA and preranked GSEA summaries for one DE result."""
    if not isinstance(result, Mapping) or result.get("available") is False:
        return None
    rows = _result_gene_table(result)
    if len(rows) < max(2, int(min_overlap)):
        return None
    universe_genes = [row["gene"] for row in rows]
    matched_sets = _match_gene_sets(
        gene_sets,
        universe_genes,
        min_overlap=int(min_overlap),
        max_size=int(max_pathway_size),
    )
    if not matched_sets:
        return None
    padj_cutoff = _finite_float(result.get("padj_cutoff")) or 0.05
    log2fc_cutoff = _finite_float(result.get("log2fc_cutoff")) or 0.5
    min_pct_expressed = _finite_float(result.get("min_pct_expressed")) or 0.0
    if min_pct_expressed > 1.0:
        min_pct_expressed /= 100.0
    ora = {
        "up": _ora_rows(
            rows,
            matched_sets,
            direction="up",
            padj_cutoff=padj_cutoff,
            log2fc_cutoff=log2fc_cutoff,
            min_pct_expressed=min_pct_expressed,
            min_overlap=int(min_overlap),
            top_n=int(top_n),
        ),
        "down": _ora_rows(
            rows,
            matched_sets,
            direction="down",
            padj_cutoff=padj_cutoff,
            log2fc_cutoff=log2fc_cutoff,
            min_pct_expressed=min_pct_expressed,
            min_overlap=int(min_overlap),
            top_n=int(top_n),
        ),
    }
    gsea = _gsea_rows(
        rows,
        matched_sets,
        min_overlap=int(min_overlap),
        top_n=int(top_n),
        permutations=int(gsea_permutations),
        seed=int(seed),
    )
    if not any(ora.values()) and not any(gsea.values()):
        return None
    return {
        "ora": ora,
        "gsea": gsea,
        "tested_pathways": int(len(matched_sets)),
        "universe_size": int(len(universe_genes)),
    }


def add_pathway_enrichment_to_pseudobulk_de(
    pseudobulk_de: Any,
    *,
    pathway_gmt: Optional[Union[str, Sequence[str]]] = None,
    top_n: int = 20,
    min_overlap: int = 3,
    max_pathway_size: int = 500,
    gsea_permutations: int = 100,
    seed: int = 0,
    organism: str = "Human",
    n_cpus: int = 1,
    progress_callback: Optional[Callable[[Mapping[str, Any]], None]] = None,
) -> Dict[str, Any]:
    """Attach pathway enrichment summaries to every available pseudobulk DE result."""
    settings: Dict[str, Any] = {
        "available": False,
        "top_n": int(top_n),
        "min_overlap": int(min_overlap),
        "max_pathway_size": int(max_pathway_size),
        "gsea_permutations": int(gsea_permutations),
        "organism": organism,
        "n_cpus": max(1, int(n_cpus)),
        "comparisons": 0,
        "enriched_comparisons": 0,
    }
    if not isinstance(pseudobulk_de, dict) or not pseudobulk_de:
        settings["reason"] = "no_pseudobulk_de"
        return settings

    comparison_items: List[Tuple[str, str, str, Dict[str, Any]]] = []
    for color_col, by_source in pseudobulk_de.items():
        if not isinstance(by_source, dict):
            continue
        for source, by_reference in by_source.items():
            if str(source).startswith("_") or not isinstance(by_reference, dict):
                continue
            for reference, result in by_reference.items():
                if str(reference).startswith("_") or not isinstance(result, dict):
                    continue
                comparison_items.append((str(color_col), str(source), str(reference), result))
    if not comparison_items:
        settings["reason"] = "no_pseudobulk_comparisons"
        return settings

    try:
        gene_sets, source_info = load_pathway_gene_sets(pathway_gmt, organism=organism)
    except Exception as exc:
        settings["reason"] = "pathway_gene_sets_unavailable"
        settings["error"] = str(exc)
        return settings
    settings.update(source_info)
    settings["available"] = True

    def _compute_comparison(item: Tuple[int, str, str, str, Dict[str, Any]]) -> Dict[str, Any]:
        submit_index, color_col, source, reference, result = item
        enrichment = enrich_pseudobulk_result(
            result,
            gene_sets,
            top_n=int(top_n),
            min_overlap=int(min_overlap),
            max_pathway_size=int(max_pathway_size),
            gsea_permutations=int(gsea_permutations),
            seed=int(seed) + submit_index,
        )
        ora = (enrichment or {}).get("ora") if isinstance(enrichment, dict) else {}
        gsea = (enrichment or {}).get("gsea") if isinstance(enrichment, dict) else {}
        return {
            "submit_index": int(submit_index),
            "color_col": color_col,
            "source": source,
            "reference": reference,
            "result": result,
            "enrichment": enrichment,
            "ora_up": len((ora or {}).get("up") or []),
            "ora_down": len((ora or {}).get("down") or []),
            "gsea_positive": len((gsea or {}).get("positive") or []),
            "gsea_negative": len((gsea or {}).get("negative") or []),
            "stored": bool(enrichment),
        }

    def _store_completed(completed: Mapping[str, Any], completion_index: int) -> bool:
        result = completed["result"]
        enrichment = completed["enrichment"]
        if enrichment:
            result["pathway_enrichment"] = enrichment
        if progress_callback is not None:
            progress_callback({
                "event": "comparison_done",
                "index": int(completion_index),
                "total": int(total_comparisons),
                "submit_index": int(completed["submit_index"]),
                "color_col": completed["color_col"],
                "source": completed["source"],
                "reference": completed["reference"],
                "ora_up": int(completed["ora_up"]),
                "ora_down": int(completed["ora_down"]),
                "gsea_positive": int(completed["gsea_positive"]),
                "gsea_negative": int(completed["gsea_negative"]),
                "stored": bool(completed["stored"]),
            })
        return bool(enrichment)

    comparison_count = 0
    enriched_count = 0
    total_comparisons = len(comparison_items)
    indexed_items = [
        (idx, color_col, source, reference, result)
        for idx, (color_col, source, reference, result) in enumerate(comparison_items, start=1)
    ]
    max_workers = min(max(1, int(n_cpus)), max(1, total_comparisons))
    if max_workers <= 1 or total_comparisons <= 1:
        for item in indexed_items:
            comparison_count += 1
            if _store_completed(_compute_comparison(item), comparison_count):
                enriched_count += 1
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_compute_comparison, item) for item in indexed_items]
            for future in as_completed(futures):
                comparison_count += 1
                if _store_completed(future.result(), comparison_count):
                    enriched_count += 1
    settings["comparisons"] = int(comparison_count)
    settings["enriched_comparisons"] = int(enriched_count)
    return settings
