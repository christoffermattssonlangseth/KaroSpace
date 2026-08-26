"""
Desktop GUI for exporting KaroSpace HTML viewers.
"""

from __future__ import annotations

import argparse
import json
import re
import threading
import traceback
from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
except Exception:  # pragma: no cover - environment dependent
    tk = None
    filedialog = None
    messagebox = None
    ttk = None

_TTK_FRAME_BASE = ttk.Frame if ttk is not None else object


def _parse_tokens(value: str) -> Optional[List[str]]:
    values = [token.strip() for token in re.split(r"[,\n;]+", str(value)) if token.strip()]
    return values or None


def _unique(values: Optional[List[str]]) -> Optional[List[str]]:
    if not values:
        return None
    seen = set()
    out: List[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out or None


def _parse_spot_size(token: str) -> Union[str, float]:
    cleaned = str(token).strip()
    if cleaned.lower() in {"auto", "adaptive", "density"}:
        return "auto"
    value = float(cleaned)
    if value <= 0:
        raise ValueError("Spot size must be a positive number or 'auto'.")
    return value


def _parse_neighbor_permutations(token: str) -> Optional[int]:
    cleaned = str(token).strip()
    if cleaned.lower() == "auto":
        return None
    value = int(cleaned)
    if value < 0:
        raise ValueError("Neighbor permutations must be >= 0 or 'auto'.")
    return value


def _parse_positive_int(name: str, token: str) -> int:
    value = int(str(token).strip())
    if value <= 0:
        raise ValueError(f"{name} must be > 0.")
    return value


def _parse_non_negative_int(name: str, token: str) -> int:
    value = int(str(token).strip())
    if value < 0:
        raise ValueError(f"{name} must be >= 0.")
    return value


FIELD_HELP_TEXT = """Input / Output
- Input (.h5ad): Path to an AnnData file. Use Browse or type a full path.
- Output (.html): Path for the exported viewer HTML file.

Core Options
- Group by: Column in adata.obs that identifies sections (for example sample_id).
- Main cells annotation: Starting cell annotation variable (obs column or a pre-loaded gene).
- Outline by: Optional obs/metadata column used to paint panel outlines. Leave blank to disable.
- Theme: light or dark.
- Title: Viewer page title shown in header/browser tab.
- Min panel size: Integer > 0 (typical 100-200).
- Spot size: "auto" or a positive number (for example 1.2).
- Downsample: Empty (disabled) or integer > 0 (max cells per section).

Dataset Loading Options
- Metadata section: Comma/newline-separated obs columns used as section metadata/filter chips.
  Empty uses defaults (course, region, condition, timepoint, ...).
- Metadata max columns: Empty or integer >= 0.
- Metadata value order (JSON): Dict mapping metadata column -> ordered values.
  Example:
  {"condition": ["naive", "treated"], "stage": ["early", "late"]}
- Metadata labels (JSON): Dict mapping column key -> display name in the viewer.
  Example:
  {"sample_id": "Sample", "last_score": "Disease score"}

Annotation & Gene Content
- Cells annotations: Comma/newline-separated obs columns to include in the annotation dropdown.
- Genes: Hand-picked gene symbols (comma/newline-separated) to pre-load for expression view.
- In embedded mode, significant shared-fit pseudobulk DE genes are embedded automatically up to the per-comparison cap.

Advanced Options
- Gene encoding: auto | dense | sparse.
- Feature storage: embedded keeps requested/top DE genes in the HTML; sidecar writes all gene expression vectors to separate sidecar files.
- Feature manifest path: Optional override for the sidecar JSON path. Leave blank to place it next to the HTML.
- Feature sparse threshold: Number between 0 and 1 (used when encoding=auto).
- Neighbor permutations: "auto" or integer >= 0.
- Neighbor stats annotations: "auto" or comma/newline-separated obs columns.
- Neighbor stats seed: Integer random seed.
- Interaction top targets: Integer > 0.
- Interaction top genes: Integer > 0.
- Interaction min cells: Integer > 0.
- Interaction min neighbors: Integer > 0.
- Pseudobulk DE uses one shared replicate + annotation fit for each selected annotation, with category-versus-category and balanced-rest contrasts. Non-embedded DE genes remain visible but cannot be clicked for expression.
- Interaction markers use the same pseudobulk counts layer and section replicate.
"""

UI_COLORS = {
    "page_bg": "#edf2f7",
    "card_bg": "#f8fafc",
    "header_bg": "#0f172a",
    "header_fg": "#f8fafc",
    "text": "#0f172a",
    "muted": "#475569",
    "border": "#cbd5e1",
    "input_bg": "#ffffff",
    "primary": "#0ea5e9",
    "primary_hover": "#0284c7",
    "primary_pressed": "#0369a1",
    "secondary": "#e2e8f0",
    "secondary_hover": "#cbd5e1",
    "secondary_pressed": "#94a3b8",
    "log_bg": "#0b1220",
    "log_fg": "#dbeafe",
}


class SearchableListEditor(_TTK_FRAME_BASE):
    def __init__(
        self,
        parent,
        title: str,
        *,
        height: int = 6,
        help_text: Optional[str] = None,
    ):
        super().__init__(parent, style="Card.TFrame")
        self._choices: List[str] = []
        self._var = tk.StringVar(value="")

        self.columnconfigure(0, weight=1)
        if title:
            ttk.Label(self, text=title).grid(row=0, column=0, sticky="w", pady=(0, 4))

        controls = ttk.Frame(self, style="Card.TFrame")
        controls.grid(row=1, column=0, sticky="ew")
        controls.columnconfigure(0, weight=1)

        self.combo = ttk.Combobox(controls, textvariable=self._var, state="normal")
        self.combo.grid(row=0, column=0, sticky="ew", padx=(0, 6))
        self.combo.bind("<KeyRelease>", self._on_type)
        self.combo.bind("<Return>", lambda _event: self.add_current())

        self.add_btn = ttk.Button(controls, text="+ Add", style="Secondary.TButton", command=self.add_current)
        self.add_btn.grid(row=0, column=1, padx=(0, 6))
        self.remove_btn = ttk.Button(controls, text="Remove", style="Secondary.TButton", command=self.remove_selected)
        self.remove_btn.grid(row=0, column=2, padx=(0, 6))
        self.clear_btn = ttk.Button(controls, text="Clear", style="Secondary.TButton", command=self.clear)
        self.clear_btn.grid(row=0, column=3)

        list_frame = ttk.Frame(self, style="Card.TFrame")
        list_frame.grid(row=2, column=0, sticky="ew", pady=(8, 0))
        list_frame.columnconfigure(0, weight=1)

        self.listbox = tk.Listbox(
            list_frame,
            height=height,
            selectmode="extended",
            activestyle="none",
            relief="solid",
            bd=1,
            highlightthickness=0,
            bg=UI_COLORS["input_bg"],
            fg=UI_COLORS["text"],
            selectbackground=UI_COLORS["primary"],
            selectforeground="#ffffff",
            font=("Helvetica", 10),
        )
        self.listbox.grid(row=0, column=0, sticky="ew")
        self.scrollbar = ttk.Scrollbar(list_frame, orient="vertical", command=self.listbox.yview)
        self.scrollbar.grid(row=0, column=1, sticky="ns")
        self.listbox.configure(yscrollcommand=self.scrollbar.set)

        if help_text:
            ttk.Label(self, text=help_text, style="Muted.TLabel").grid(row=3, column=0, sticky="w", pady=(6, 0))

    def _on_type(self, _event) -> None:
        self._update_combo_values(self._var.get().strip())

    def _update_combo_values(self, query: str = "") -> None:
        q = str(query or "").strip().lower()
        if not q:
            filtered = self._choices
        else:
            filtered = [v for v in self._choices if q in v.lower()]
        self.combo.configure(values=filtered[:300])

    def set_choices(self, values: Sequence[str]) -> None:
        self._choices = sorted({str(v) for v in values if str(v).strip()})
        self._update_combo_values(self._var.get())

    def add_item(self, value: str) -> None:
        token = str(value).strip()
        if not token:
            return
        existing = self.get_items()
        if token in existing:
            idx = existing.index(token)
            self.listbox.selection_clear(0, "end")
            self.listbox.selection_set(idx)
            self.listbox.see(idx)
            return
        self.listbox.insert("end", token)

    def add_current(self) -> None:
        token = self._var.get().strip()
        if not token:
            return
        self.add_item(token)
        self._var.set("")
        self._update_combo_values("")

    def remove_selected(self) -> None:
        indices = list(self.listbox.curselection())
        for idx in reversed(indices):
            self.listbox.delete(idx)

    def clear(self) -> None:
        self.listbox.delete(0, "end")

    def get_items(self) -> List[str]:
        return [str(item) for item in self.listbox.get(0, "end")]

    def set_items(self, values: Sequence[str]) -> None:
        self.clear()
        for value in _unique([str(v).strip() for v in values if str(v).strip()]) or []:
            self.listbox.insert("end", value)

    def set_enabled(self, enabled: bool) -> None:
        state = "normal" if enabled else "disabled"
        self.combo.configure(state=state)
        self.add_btn.configure(state=state)
        self.remove_btn.configure(state=state)
        self.clear_btn.configure(state=state)
        self.listbox.configure(state=state)


class KaroSpaceExportGUI:
    def __init__(self, root: tk.Tk, initial_input: Optional[str] = None, initial_output: Optional[str] = None):
        self.root = root
        self.root.title("KaroSpaceBuilder")
        self.root.geometry("1120x920")
        self.root.minsize(1000, 760)
        self._configure_styles()

        self._obs_columns: List[str] = []
        self._var_names: List[str] = []
        self._obsm_keys: List[str] = []
        self._is_busy = False

        self.input_path = tk.StringVar(value=initial_input or "")
        self.output_path = tk.StringVar(value=initial_output or "karospace.html")
        self.section_key = tk.StringVar(value="sample_id")
        self.spatial_key = tk.StringVar(value="spatial")
        self.annotation = tk.StringVar(value="leiden")
        self.outline_by = tk.StringVar(value="condition")
        self.title = tk.StringVar(value="KaroSpace")
        self.section_metadata = tk.StringVar(value="condition")
        self.metadata_max_columns = tk.StringVar(value="")
        self.metadata_labels = tk.StringVar(value='{\n  "condition": "Condition"\n}')
        self.min_panel_size = tk.StringVar(value="150")
        self.spot_size = tk.StringVar(value="auto")
        self.downsample = tk.StringVar(value="")
        self.gene_correlation_top_n = tk.StringVar(value="5")
        self.spatial_variable_genes_n = tk.StringVar(value="200")
        self.feature_encoding = tk.StringVar(value="auto")
        self.feature_storage = tk.StringVar(value="embedded")
        self.feature_manifest_path = tk.StringVar(value="")
        self.feature_sparse_zero_threshold = tk.StringVar(value="0.8")
        self.neighbor_permutations = tk.StringVar(value="auto")
        self.neighbor_stats_seed = tk.StringVar(value="0")
        self.interaction_markers_top_targets = tk.StringVar(value="5")
        self.interaction_markers_top_genes = tk.StringVar(value="20")
        self.interaction_markers_min_cells = tk.StringVar(value="30")
        self.interaction_markers_min_neighbors = tk.StringVar(value="1")
        self.status_text = tk.StringVar(value="Ready.")
        self.advanced_expanded = tk.BooleanVar(value=False)
        self.neighbor_stats_auto = tk.BooleanVar(value=True)

        self._build_layout()
        self._bind_shortcuts()

        if initial_input:
            self._prefill_output_from_input(initial_input)
            self.inspect_dataset()

    def _configure_styles(self) -> None:
        self.root.configure(bg=UI_COLORS["page_bg"])
        self.root.option_add("*tearOff", False)

        style = ttk.Style(self.root)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass

        style.configure(
            ".",
            background=UI_COLORS["page_bg"],
            foreground=UI_COLORS["text"],
            font=("Helvetica", 11),
        )
        style.configure("App.TFrame", background=UI_COLORS["page_bg"])
        style.configure("Card.TFrame", background=UI_COLORS["card_bg"])
        style.configure("Header.TFrame", background=UI_COLORS["header_bg"])
        style.configure(
            "AppTitle.TLabel",
            background=UI_COLORS["header_bg"],
            foreground=UI_COLORS["header_fg"],
            font=("Helvetica", 24, "bold"),
        )
        style.configure(
            "Subtitle.TLabel",
            background=UI_COLORS["header_bg"],
            foreground="#cbd5e1",
            font=("Helvetica", 11),
        )
        style.configure(
            "Muted.TLabel",
            background=UI_COLORS["card_bg"],
            foreground=UI_COLORS["muted"],
            font=("Helvetica", 10),
        )
        style.configure(
            "Status.TLabel",
            background=UI_COLORS["card_bg"],
            foreground=UI_COLORS["muted"],
            font=("Helvetica", 10, "italic"),
        )
        style.configure(
            "TLabel",
            background=UI_COLORS["card_bg"],
            foreground=UI_COLORS["text"],
            font=("Helvetica", 10),
        )
        style.configure(
            "Card.TLabelframe",
            background=UI_COLORS["card_bg"],
            bordercolor=UI_COLORS["border"],
            borderwidth=1,
            relief="solid",
        )
        style.configure(
            "Card.TLabelframe.Label",
            background=UI_COLORS["card_bg"],
            foreground=UI_COLORS["text"],
            font=("Helvetica", 11, "bold"),
        )
        style.configure(
            "TEntry",
            fieldbackground=UI_COLORS["input_bg"],
            background=UI_COLORS["input_bg"],
            foreground=UI_COLORS["text"],
            bordercolor=UI_COLORS["border"],
            lightcolor=UI_COLORS["border"],
            darkcolor=UI_COLORS["border"],
            insertcolor=UI_COLORS["text"],
            padding=(8, 5),
        )
        style.configure(
            "TCombobox",
            fieldbackground=UI_COLORS["input_bg"],
            background=UI_COLORS["input_bg"],
            foreground=UI_COLORS["text"],
            bordercolor=UI_COLORS["border"],
            lightcolor=UI_COLORS["border"],
            darkcolor=UI_COLORS["border"],
            padding=(8, 5),
        )
        style.map(
            "TCombobox",
            fieldbackground=[("readonly", UI_COLORS["input_bg"])],
            selectbackground=[("readonly", UI_COLORS["input_bg"])],
            selectforeground=[("readonly", UI_COLORS["text"])],
        )
        style.configure(
            "Primary.TButton",
            background=UI_COLORS["primary"],
            foreground="#ffffff",
            borderwidth=0,
            focusthickness=0,
            padding=(12, 8),
            font=("Helvetica", 10, "bold"),
        )
        style.map(
            "Primary.TButton",
            background=[
                ("disabled", "#94a3b8"),
                ("pressed", UI_COLORS["primary_pressed"]),
                ("active", UI_COLORS["primary_hover"]),
            ],
            foreground=[("disabled", "#e2e8f0")],
        )
        style.configure(
            "Secondary.TButton",
            background=UI_COLORS["secondary"],
            foreground=UI_COLORS["text"],
            borderwidth=0,
            focusthickness=0,
            padding=(12, 8),
            font=("Helvetica", 10, "bold"),
        )
        style.map(
            "Secondary.TButton",
            background=[
                ("disabled", "#e2e8f0"),
                ("pressed", UI_COLORS["secondary_pressed"]),
                ("active", UI_COLORS["secondary_hover"]),
            ],
        )
        style.configure(
            "TCheckbutton",
            background=UI_COLORS["card_bg"],
            foreground=UI_COLORS["text"],
            font=("Helvetica", 10),
        )
        style.configure("TNotebook", background=UI_COLORS["page_bg"], borderwidth=0)
        style.configure(
            "TNotebook.Tab",
            background=UI_COLORS["secondary"],
            foreground=UI_COLORS["muted"],
            padding=(14, 8),
            font=("Helvetica", 10, "bold"),
        )
        style.map(
            "TNotebook.Tab",
            background=[
                ("selected", UI_COLORS["card_bg"]),
                ("active", UI_COLORS["secondary_hover"]),
            ],
            foreground=[
                ("selected", UI_COLORS["text"]),
                ("active", UI_COLORS["text"]),
            ],
        )
        style.configure("TScrollbar", background=UI_COLORS["secondary"], troughcolor=UI_COLORS["page_bg"])

    def _build_layout(self) -> None:
        shell = ttk.Frame(self.root, style="App.TFrame")
        shell.pack(fill="both", expand=True)

        self.main_canvas = tk.Canvas(
            shell,
            highlightthickness=0,
            borderwidth=0,
            background=UI_COLORS["page_bg"],
        )
        self.main_scrollbar = ttk.Scrollbar(shell, orient="vertical", command=self.main_canvas.yview)
        self.main_canvas.configure(yscrollcommand=self.main_scrollbar.set)

        self.main_scrollbar.pack(side="right", fill="y")
        self.main_canvas.pack(side="left", fill="both", expand=True)

        container = ttk.Frame(self.main_canvas, style="App.TFrame", padding=16)
        self._canvas_window = self.main_canvas.create_window((0, 0), window=container, anchor="nw")
        container.bind("<Configure>", self._on_main_content_configure)
        self.main_canvas.bind("<Configure>", self._on_main_canvas_configure)
        self.main_canvas.bind_all("<MouseWheel>", self._on_mousewheel)
        self.main_canvas.bind_all("<Button-4>", self._on_mousewheel)
        self.main_canvas.bind_all("<Button-5>", self._on_mousewheel)

        header = ttk.Frame(container, style="Header.TFrame", padding=(18, 14))
        header.pack(fill="x", pady=(0, 12))
        ttk.Label(header, text="KaroSpaceBuilder", style="AppTitle.TLabel").pack(anchor="w")
        ttk.Label(
            header,
            text="Build and export shareable KaroSpace viewers with full control over annotations, genes, and analytics.",
            style="Subtitle.TLabel",
        ).pack(anchor="w", pady=(4, 0))

        preset_group = ttk.LabelFrame(container, text="Presets", padding=10, style="Card.TLabelframe")
        preset_group.pack(fill="x", pady=(0, 10))
        preset_row = ttk.Frame(preset_group, style="Card.TFrame")
        preset_row.pack(fill="x")
        ttk.Button(
            preset_row,
            text="Default",
            style="Secondary.TButton",
            command=lambda: self._apply_preset("default"),
        ).pack(side="left")
        ttk.Button(
            preset_row,
            text="Pancreas",
            style="Secondary.TButton",
            command=lambda: self._apply_preset("pancreas"),
        ).pack(side="left", padx=(8, 0))
        ttk.Button(
            preset_row,
            text="Lightweight",
            style="Secondary.TButton",
            command=lambda: self._apply_preset("lightweight"),
        ).pack(side="left", padx=(8, 0))
        ttk.Label(
            preset_row,
            text="Apply a profile to prefill fields instantly.",
            style="Muted.TLabel",
        ).pack(side="left", padx=(12, 0))

        self.tabs = ttk.Notebook(container)
        self.tabs.pack(fill="both", expand=True, pady=(0, 10))

        basic_tab = ttk.Frame(self.tabs, style="App.TFrame", padding=4)
        colors_tab = ttk.Frame(self.tabs, style="App.TFrame", padding=4)
        advanced_tab = ttk.Frame(self.tabs, style="App.TFrame", padding=4)
        help_tab = ttk.Frame(self.tabs, style="App.TFrame", padding=4)
        self.tabs.add(basic_tab, text="Basic")
        self.tabs.add(colors_tab, text="Annotations & Genes")
        self.tabs.add(advanced_tab, text="Advanced")
        self.tabs.add(help_tab, text="Help")

        io_group = ttk.LabelFrame(basic_tab, text="Input / Output", padding=12, style="Card.TLabelframe")
        io_group.pack(fill="x", pady=(0, 10))
        io_group.columnconfigure(1, weight=1)
        ttk.Label(io_group, text="Input (.h5ad)").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=4)
        ttk.Entry(io_group, textvariable=self.input_path).grid(row=0, column=1, sticky="ew", pady=4)
        ttk.Button(io_group, text="Browse", command=self.choose_input, style="Secondary.TButton").grid(
            row=0, column=2, padx=(8, 0), pady=4
        )
        ttk.Button(io_group, text="Inspect", command=self.inspect_dataset, style="Secondary.TButton").grid(
            row=0, column=3, padx=(8, 0), pady=4
        )
        ttk.Label(io_group, text="Output (.html)").grid(row=1, column=0, sticky="w", padx=(0, 8), pady=4)
        ttk.Entry(io_group, textvariable=self.output_path).grid(row=1, column=1, sticky="ew", pady=4)
        ttk.Button(io_group, text="Browse", command=self.choose_output, style="Secondary.TButton").grid(
            row=1, column=2, padx=(8, 0), pady=4
        )

        core_group = ttk.LabelFrame(basic_tab, text="Core Options", padding=12, style="Card.TLabelframe")
        core_group.pack(fill="x", pady=(0, 10))
        core_group.columnconfigure(1, weight=1)
        core_group.columnconfigure(3, weight=1)
        ttk.Label(core_group, text="Section key").grid(row=0, column=0, sticky="w", pady=4)
        self.section_key_combo = ttk.Combobox(core_group, textvariable=self.section_key, state="normal")
        self.section_key_combo.grid(row=0, column=1, sticky="ew", padx=(8, 16), pady=4)
        ttk.Label(core_group, text="Main cells annotation").grid(row=0, column=2, sticky="w", pady=4)
        self.color_combo = ttk.Combobox(core_group, textvariable=self.annotation, state="normal")
        self.color_combo.grid(row=0, column=3, sticky="ew", pady=4)
        ttk.Label(core_group, text="Outline by").grid(row=1, column=0, sticky="w", pady=4)
        self.outline_combo = ttk.Combobox(core_group, textvariable=self.outline_by, state="normal")
        self.outline_combo.grid(row=1, column=1, sticky="ew", padx=(8, 16), pady=4)
        ttk.Label(core_group, text="Title").grid(row=1, column=2, sticky="w", pady=4)
        ttk.Entry(core_group, textvariable=self.title).grid(row=1, column=3, sticky="ew", pady=4)
        ttk.Label(core_group, text="Min panel size").grid(row=2, column=0, sticky="w", pady=4)
        ttk.Entry(core_group, textvariable=self.min_panel_size).grid(row=2, column=1, sticky="ew", padx=(8, 16), pady=4)
        ttk.Label(core_group, text="Spot size").grid(row=2, column=2, sticky="w", pady=4)
        ttk.Entry(core_group, textvariable=self.spot_size).grid(row=2, column=3, sticky="ew", pady=4)
        ttk.Label(core_group, text="Downsample").grid(row=3, column=0, sticky="w", pady=4)
        ttk.Entry(core_group, textvariable=self.downsample).grid(row=3, column=1, sticky="ew", padx=(8, 16), pady=4)

        dataset_group = ttk.LabelFrame(basic_tab, text="Dataset Loading Options", padding=12, style="Card.TLabelframe")
        dataset_group.pack(fill="x", pady=(0, 10))
        dataset_group.columnconfigure(1, weight=1)
        dataset_group.columnconfigure(3, weight=1)
        ttk.Label(dataset_group, text="Metadata section").grid(row=0, column=0, sticky="w", pady=4)
        ttk.Entry(dataset_group, textvariable=self.section_metadata).grid(row=0, column=1, sticky="ew", padx=(8, 16), pady=4)
        ttk.Label(dataset_group, text="Metadata max columns").grid(row=0, column=2, sticky="w", pady=4)
        ttk.Entry(dataset_group, textvariable=self.metadata_max_columns).grid(row=0, column=3, sticky="ew", pady=4)
        ttk.Label(dataset_group, text="Metadata value order (JSON)").grid(row=1, column=0, sticky="nw", pady=4)
        self.metadata_value_order_text = tk.Text(dataset_group, height=4, wrap="word")
        self.metadata_value_order_text.grid(row=1, column=1, columnspan=3, sticky="ew", pady=4)
        self.metadata_value_order_text.insert("1.0", '{\n  "condition": ["control", "treated"]\n}')
        self.metadata_value_order_text.configure(
            bg=UI_COLORS["input_bg"],
            fg=UI_COLORS["text"],
            insertbackground=UI_COLORS["text"],
            relief="solid",
            bd=1,
            highlightthickness=0,
            padx=8,
            pady=6,
            font=("Helvetica", 10),
        )

        ttk.Label(dataset_group, text="Metadata labels (JSON)").grid(row=2, column=0, sticky="nw", pady=4)
        self.metadata_labels_text = tk.Text(dataset_group, height=4, wrap="word")
        self.metadata_labels_text.grid(row=2, column=1, columnspan=3, sticky="ew", pady=4)
        self.metadata_labels_text.insert("1.0", self.metadata_labels.get())
        self.metadata_labels_text.configure(
            bg=UI_COLORS["input_bg"],
            fg=UI_COLORS["text"],
            insertbackground=UI_COLORS["text"],
            relief="solid",
            bd=1,
            highlightthickness=0,
            padx=8,
            pady=6,
            font=("Helvetica", 10),
        )

        ttk.Label(dataset_group, text="Spatial key").grid(row=3, column=0, sticky="w", pady=4)
        self.spatial_key_combo = ttk.Combobox(dataset_group, textvariable=self.spatial_key, state="normal")
        self.spatial_key_combo.grid(row=3, column=1, sticky="ew", padx=(8, 16), pady=4)


        colors_group = ttk.LabelFrame(colors_tab, text="Annotation Layers", padding=12, style="Card.TLabelframe")
        colors_group.pack(fill="x", pady=(0, 10))
        self.cell_annotations_editor = SearchableListEditor(
            colors_group,
            "Cells annotations",
            height=6,
            help_text="Search obs columns, then click + Add. This controls extra items in the viewer annotation dropdown.",
        )
        self.cell_annotations_editor.pack(fill="x")

        genes_group = ttk.LabelFrame(colors_tab, text="Gene Layers", padding=12, style="Card.TLabelframe")
        genes_group.pack(fill="x", pady=(0, 10))
        self.genes_editor = SearchableListEditor(
            genes_group,
            "Genes",
            height=7,
            help_text="Hand-pick genes for expression view and dotplot.",
        )
        self.genes_editor.pack(fill="x")

        genes_opts = ttk.Frame(genes_group, style="Card.TFrame")
        genes_opts.pack(fill="x", pady=(10, 0))
        ttk.Label(genes_opts, text="Corr. top N").pack(side="left", padx=(0, 4))
        ttk.Entry(genes_opts, textvariable=self.gene_correlation_top_n, width=6).pack(side="left")
        ttk.Label(genes_opts, text="Moran genes N").pack(side="left", padx=(12, 4))
        ttk.Entry(genes_opts, textvariable=self.spatial_variable_genes_n, width=6).pack(side="left")

        advanced_group = ttk.LabelFrame(advanced_tab, text="Advanced Options", padding=12, style="Card.TLabelframe")
        advanced_group.pack(fill="x", pady=(0, 10))
        advanced_group.columnconfigure(0, weight=1)
        advanced_header = ttk.Frame(advanced_group, style="Card.TFrame")
        advanced_header.grid(row=0, column=0, sticky="ew")
        self.advanced_toggle_btn = ttk.Button(
            advanced_header,
            text="Show Advanced Options",
            style="Secondary.TButton",
            command=self._toggle_advanced_options,
        )
        self.advanced_toggle_btn.pack(anchor="w")
        ttk.Label(
            advanced_header,
            text="Expand to tune encoding, neighbor stats, and interaction marker settings.",
            style="Muted.TLabel",
        ).pack(anchor="w", pady=(6, 0))

        self.advanced_content = ttk.Frame(advanced_group, style="Card.TFrame")
        self.advanced_content.grid(row=1, column=0, sticky="ew", pady=(10, 0))

        encoding_group = ttk.LabelFrame(self.advanced_content, text="Encoding", padding=10, style="Card.TLabelframe")
        encoding_group.pack(fill="x", pady=(0, 10))
        encoding_group.columnconfigure(1, weight=1)
        encoding_group.columnconfigure(3, weight=1)
        ttk.Label(encoding_group, text="Gene encoding").grid(row=0, column=0, sticky="w", pady=4)
        ttk.Combobox(encoding_group, textvariable=self.feature_encoding, values=["auto", "dense", "sparse"], state="readonly").grid(
            row=0, column=1, sticky="ew", padx=(8, 16), pady=4
        )
        ttk.Label(encoding_group, text="Feature storage").grid(row=0, column=2, sticky="w", pady=4)
        ttk.Combobox(encoding_group, textvariable=self.feature_storage, values=["embedded", "sidecar"], state="readonly").grid(
            row=0, column=3, sticky="ew", pady=4
        )
        ttk.Label(encoding_group, text="Feature sparse threshold").grid(row=1, column=0, sticky="w", pady=4)
        ttk.Entry(encoding_group, textvariable=self.feature_sparse_zero_threshold).grid(row=1, column=1, sticky="ew", padx=(8, 16), pady=4)
        ttk.Label(encoding_group, text="Feature manifest path").grid(row=1, column=2, sticky="w", pady=4)
        ttk.Entry(encoding_group, textvariable=self.feature_manifest_path).grid(row=1, column=3, sticky="ew", pady=4)

        neighbor_group = ttk.LabelFrame(self.advanced_content, text="Neighbor Stats", padding=10, style="Card.TLabelframe")
        neighbor_group.pack(fill="x", pady=(0, 10))
        neighbor_group.columnconfigure(1, weight=1)
        neighbor_group.columnconfigure(3, weight=1)
        ttk.Label(neighbor_group, text="Neighbor permutations").grid(row=0, column=0, sticky="w", pady=4)
        ttk.Entry(neighbor_group, textvariable=self.neighbor_permutations).grid(row=0, column=1, sticky="ew", padx=(8, 16), pady=4)
        ttk.Label(neighbor_group, text="Neighbor stats seed").grid(row=0, column=2, sticky="w", pady=4)
        ttk.Entry(neighbor_group, textvariable=self.neighbor_stats_seed).grid(row=0, column=3, sticky="ew", pady=4)
        ttk.Checkbutton(
            neighbor_group,
            text="Auto neighbor annotations (use main cell annotation)",
            variable=self.neighbor_stats_auto,
            command=self._on_neighbor_auto_toggle,
        ).grid(row=1, column=0, columnspan=4, sticky="w", pady=(2, 8))
        self.neighbor_stats_annotations_editor = SearchableListEditor(
            neighbor_group,
            "Neighbor stats annotation columns",
            height=5,
        )
        self.neighbor_stats_annotations_editor.grid(row=2, column=0, columnspan=4, sticky="ew")

        interaction_group = ttk.LabelFrame(self.advanced_content, text="Interaction Markers", padding=10, style="Card.TLabelframe")
        interaction_group.pack(fill="x")
        interaction_group.columnconfigure(1, weight=1)
        interaction_group.columnconfigure(3, weight=1)
        ttk.Label(interaction_group, text="Top targets").grid(row=0, column=0, sticky="w", pady=4)
        ttk.Entry(interaction_group, textvariable=self.interaction_markers_top_targets).grid(
            row=0, column=1, sticky="ew", padx=(8, 16), pady=4
        )
        ttk.Label(interaction_group, text="Top genes").grid(row=0, column=2, sticky="w", pady=4)
        ttk.Entry(interaction_group, textvariable=self.interaction_markers_top_genes).grid(row=0, column=3, sticky="ew", pady=4)
        ttk.Label(interaction_group, text="Min cells").grid(row=1, column=0, sticky="w", pady=4)
        ttk.Entry(interaction_group, textvariable=self.interaction_markers_min_cells).grid(
            row=1, column=1, sticky="ew", padx=(8, 16), pady=4
        )
        ttk.Label(interaction_group, text="Min neighbors").grid(row=1, column=2, sticky="w", pady=4)
        ttk.Entry(interaction_group, textvariable=self.interaction_markers_min_neighbors).grid(
            row=1, column=3, sticky="ew", pady=4
        )
        self._set_advanced_options_visible(False)
        self._on_neighbor_auto_toggle()

        hint_group = ttk.LabelFrame(help_tab, text="Hints", padding=12, style="Card.TLabelframe")
        hint_group.pack(fill="x")
        hint_lines = [
            "1. Start with a preset (Default, Pancreas, Lightweight).",
            "2. Browse to your .h5ad and click Inspect.",
            "3. Use tabbed editors to add annotations, genes, and section metadata lists.",
            "4. Open Advanced only when needed.",
            "5. Export HTML.",
        ]
        ttk.Label(hint_group, text="\n".join(hint_lines), justify="left", style="Muted.TLabel").pack(anchor="w")

        guide_group = ttk.LabelFrame(help_tab, text="Field Guide", padding=10, style="Card.TLabelframe")
        guide_group.pack(fill="both", expand=True, pady=(10, 0))
        guide_group.rowconfigure(0, weight=1)
        guide_group.columnconfigure(0, weight=1)
        guide_group.columnconfigure(1, weight=0)
        guide_text = tk.Text(guide_group, wrap="word", height=32)
        guide_text.grid(row=0, column=0, sticky="nsew")
        guide_scroll = ttk.Scrollbar(guide_group, orient="vertical", command=guide_text.yview)
        guide_scroll.grid(row=0, column=1, sticky="ns")
        guide_text.configure(yscrollcommand=guide_scroll.set)
        guide_text.insert("1.0", FIELD_HELP_TEXT)
        guide_text.configure(
            state="disabled",
            bg=UI_COLORS["input_bg"],
            fg=UI_COLORS["text"],
            relief="solid",
            bd=1,
            highlightthickness=0,
            padx=10,
            pady=10,
            font=("Helvetica", 10),
        )

        actions = ttk.Frame(container, style="Card.TFrame", padding=(10, 10))
        actions.pack(fill="x", pady=(0, 10))
        self.inspect_btn = ttk.Button(
            actions,
            text="Inspect Dataset",
            command=self.inspect_dataset,
            style="Secondary.TButton",
        )
        self.inspect_btn.pack(side="left")
        self.export_btn = ttk.Button(
            actions,
            text="Export HTML",
            command=self.export_html,
            style="Primary.TButton",
        )
        self.export_btn.pack(side="left", padx=(8, 0))
        ttk.Label(actions, textvariable=self.status_text, style="Status.TLabel").pack(side="left", padx=(14, 0))

        log_group = ttk.LabelFrame(container, text="Log", padding=10, style="Card.TLabelframe")
        log_group.pack(fill="both", expand=True, pady=(0, 0))
        self.log_text = tk.Text(log_group, height=10, wrap="word")
        self.log_text.pack(fill="both", expand=True)
        self.log_text.configure(
            state="disabled",
            bg=UI_COLORS["log_bg"],
            fg=UI_COLORS["log_fg"],
            insertbackground=UI_COLORS["log_fg"],
            relief="solid",
            bd=1,
            highlightthickness=0,
            padx=10,
            pady=8,
            font=("Menlo", 10) if self.root.tk.call("tk", "windowingsystem") == "aqua" else ("Consolas", 10),
        )

        self._apply_preset("default", log=False)
        self._on_main_content_configure(None)

    def _toggle_advanced_options(self) -> None:
        self._set_advanced_options_visible(not bool(self.advanced_expanded.get()))

    def _set_advanced_options_visible(self, visible: bool) -> None:
        self.advanced_expanded.set(bool(visible))
        if visible:
            self.advanced_content.grid()
            self.advanced_toggle_btn.configure(text="Hide Advanced Options")
        else:
            self.advanced_content.grid_remove()
            self.advanced_toggle_btn.configure(text="Show Advanced Options")
        self._on_main_content_configure(None)

    def _on_neighbor_auto_toggle(self) -> None:
        auto_mode = bool(self.neighbor_stats_auto.get())
        if hasattr(self, "neighbor_stats_annotations_editor"):
            self.neighbor_stats_annotations_editor.set_enabled(not auto_mode)

    def _set_text_widget(self, widget: tk.Text, content: str) -> None:
        widget.delete("1.0", "end")
        widget.insert("1.0", str(content))

    @staticmethod
    def _parse_json_object_or_none(text: str):
        """Parse `text` as a JSON object, returning None if it isn't one."""
        stripped = str(text or "").strip()
        if not stripped:
            return {}
        try:
            parsed = json.loads(stripped)
        except json.JSONDecodeError:
            return None
        return parsed if isinstance(parsed, dict) else None

    def _merge_json_text_widget(self, widget: tk.Text, content: str) -> None:
        """Merge JSON object `content` into the widget's existing JSON.

        Existing keys win on conflict so manual edits are never clobbered, and
        any keys the user already added are preserved. Falls back to plain
        overwrite if either side isn't a valid JSON object.
        """
        incoming = self._parse_json_object_or_none(content)
        existing = self._parse_json_object_or_none(widget.get("1.0", "end"))
        if incoming is None or existing is None:
            self._set_text_widget(widget, content)
            return
        merged = dict(incoming)
        merged.update(existing)
        self._set_text_widget(widget, json.dumps(merged, indent=2))

    def _apply_preset(self, preset: str, *, log: bool = True) -> None:
        name = str(preset or "").strip().lower()

        common = {
            "title": "KaroSpace",
            "section_key": "sample_id",
            "spatial_key": "spatial",
            "outline_by": "condition",
            "spot_size": "auto",
            "feature_encoding": "auto",
            "feature_storage": "embedded",
            "feature_manifest_path": "",
            "feature_sparse_zero_threshold": "0.8",
        }
        for key, value in common.items():
            attr = getattr(self, key)
            if isinstance(attr, tk.BooleanVar):
                attr.set(bool(value))
            else:
                attr.set(str(value))

        if name == "pancreas":
            self.section_key.set("sample_id")
            self.annotation.set("leiden_2")
            self.outline_by.set("condition")
            self.title.set("KaroSpace")
            self.min_panel_size.set("120")
            self.downsample.set("1000000")
            self.section_metadata.set("condition")
            self.metadata_max_columns.set("")
            self._set_text_widget(self.metadata_value_order_text, '{\n  "stage": []\n}')
            self._merge_json_text_widget(self.metadata_labels_text, '{\n  "condition": "Condition",\n  "stage": "Stage"\n}')

            pancreas_colors = [
                "leiden_0.5",
                "leiden_1",
                "leiden_1.5",
                "leiden_2",
                "gmm_mana_5",
                "gmm_mana_8",
                "gmm_mana_10",
                "gmm_mana_12",
                "gmm_mana_15",
                "gmm_mana_20",
            ]
            pancreas_genes = [
                "Arg1",
                "Cd74",
                "Cldn11",
                "Col1a2",
                "Ctss",
                "Foxp3",
                "Gfap",
                "Gpnmb",
                "Grn",
                "H2-Aa",
                "H2-Ab1",
                "H2-Eb1",
                "Mbp",
                "Meg3",
                "Mki67",
                "Ptgds",
                "Serpina3n",
            ]
            self.cell_annotations_editor.set_items(pancreas_colors)
            self.genes_editor.set_items(pancreas_genes)

            self.neighbor_stats_auto.set(False)
            self.neighbor_stats_annotations_editor.set_items(["leiden_2"])
            self.neighbor_permutations.set("25")
            self.neighbor_stats_seed.set("42")

            self.interaction_markers_top_targets.set("6")
            self.interaction_markers_top_genes.set("15")
            self.interaction_markers_min_cells.set("30")
            self.interaction_markers_min_neighbors.set("1")

            self.status_text.set("Preset applied: Pancreas")
            label = "Pancreas"
        elif name == "lightweight":
            self.section_key.set("sample_id")
            self.annotation.set("leiden")
            self.outline_by.set("")
            self.title.set("KaroSpace")
            self.min_panel_size.set("140")
            self.downsample.set("50000")
            self.section_metadata.set("condition")
            self.metadata_max_columns.set("2")
            self._set_text_widget(self.metadata_value_order_text, '{\n  "condition": ["control", "treated"]\n}')
            self._merge_json_text_widget(self.metadata_labels_text, '{\n  "condition": "Condition"\n}')

            self.cell_annotations_editor.set_items(["leiden_1"])
            self.genes_editor.set_items(["Cd4", "Cd8a", "Mki67"])

            self.neighbor_stats_auto.set(True)
            self.neighbor_stats_annotations_editor.set_items([])
            self.neighbor_permutations.set("0")
            self.neighbor_stats_seed.set("0")

            self.interaction_markers_top_targets.set("4")
            self.interaction_markers_top_genes.set("10")
            self.interaction_markers_min_cells.set("20")
            self.interaction_markers_min_neighbors.set("1")

            self.status_text.set("Preset applied: Lightweight")
            label = "Lightweight"
        else:
            self.section_key.set("sample_id")
            self.annotation.set("leiden")
            self.outline_by.set("condition")
            self.title.set("KaroSpace")
            self.min_panel_size.set("150")
            self.downsample.set("")
            self.section_metadata.set("condition")
            self.metadata_max_columns.set("")
            self._set_text_widget(self.metadata_value_order_text, '{\n  "condition": ["control", "treated"]\n}')
            self._merge_json_text_widget(self.metadata_labels_text, '{\n  "condition": "Condition"\n}')

            self.cell_annotations_editor.set_items(["leiden_1", "leiden_2", "gmm_mana_10"])
            self.genes_editor.set_items(["Cd4", "Cd8a", "Gfap", "Mki67"])

            self.neighbor_stats_auto.set(True)
            self.neighbor_stats_annotations_editor.set_items([])
            self.neighbor_permutations.set("auto")
            self.neighbor_stats_seed.set("0")

            self.interaction_markers_top_targets.set("5")
            self.interaction_markers_top_genes.set("20")
            self.interaction_markers_min_cells.set("30")
            self.interaction_markers_min_neighbors.set("1")

            self.status_text.set("Preset applied: Default")
            label = "Default"

        self._on_neighbor_auto_toggle()
        self._on_main_content_configure(None)
        if log:
            self._log(f"Applied preset: {label}")

    def _on_main_content_configure(self, _event) -> None:
        if not hasattr(self, "main_canvas"):
            return
        self.main_canvas.configure(scrollregion=self.main_canvas.bbox("all"))

    def _on_main_canvas_configure(self, event) -> None:
        if not hasattr(self, "_canvas_window"):
            return
        self.main_canvas.itemconfigure(self._canvas_window, width=event.width)

    def _on_mousewheel(self, event) -> None:
        # Keep text-entry widgets' native wheel behavior unchanged.
        widget = event.widget
        if isinstance(widget, (tk.Text, tk.Entry, ttk.Entry, ttk.Combobox, tk.Listbox)):
            return

        if getattr(event, "num", None) == 4:
            delta = -1
        elif getattr(event, "num", None) == 5:
            delta = 1
        else:
            raw = int(getattr(event, "delta", 0))
            if raw == 0:
                return
            if abs(raw) >= 120:
                delta = -int(raw / 120)
            else:
                delta = -1 if raw > 0 else 1

        self.main_canvas.yview_scroll(delta, "units")

    def _bind_shortcuts(self) -> None:
        self.root.bind("<Command-Return>", lambda _: self.export_html())
        self.root.bind("<Control-Return>", lambda _: self.export_html())

    def _set_busy(self, busy: bool, status: str) -> None:
        self._is_busy = busy
        state = "disabled" if busy else "normal"
        self.inspect_btn.configure(state=state)
        self.export_btn.configure(state=state)
        self.status_text.set(status)
        self.root.configure(cursor="watch" if busy else "")
        self.root.update_idletasks()

    def _log(self, message: str) -> None:
        self.log_text.configure(state="normal")
        self.log_text.insert("end", f"{message}\n")
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    def choose_input(self) -> None:
        if self._is_busy:
            return
        assert filedialog is not None
        filename = filedialog.askopenfilename(
            title="Select input h5ad file",
            filetypes=[("AnnData files", "*.h5ad"), ("All files", "*.*")],
        )
        if not filename:
            return
        self.input_path.set(filename)
        self._prefill_output_from_input(filename)
        self.inspect_dataset()

    def choose_output(self) -> None:
        if self._is_busy:
            return
        assert filedialog is not None
        suggested = self.output_path.get().strip() or "karospace.html"
        filename = filedialog.asksaveasfilename(
            title="Choose output HTML file",
            defaultextension=".html",
            initialfile=Path(suggested).name,
            filetypes=[("HTML files", "*.html"), ("All files", "*.*")],
        )
        if filename:
            self.output_path.set(filename)

    def _prefill_output_from_input(self, input_path: str) -> None:
        current = self.output_path.get().strip()
        if current and current not in {"karospace.html", "viewer.html"}:
            return
        source = Path(input_path)
        self.output_path.set(str(source.with_suffix(".html")))

    def inspect_dataset(self) -> None:
        if self._is_busy:
            return
        assert messagebox is not None
        input_path = self.input_path.get().strip()
        if not input_path:
            messagebox.showerror("Missing input", "Choose an input .h5ad file first.")
            return
        source = Path(input_path).expanduser()
        if not source.exists():
            messagebox.showerror("Missing input", f"Input file was not found:\n{source}")
            return

        self._set_busy(True, "Inspecting dataset...")
        self._log(f"Inspecting dataset: {source}")

        def worker() -> None:
            try:
                import scanpy as sc

                adata = sc.read_h5ad(str(source), backed="r")
                obs_columns = [str(c) for c in adata.obs.columns]
                var_names = [str(v) for v in adata.var_names]
                obsm_keys = [str(k) for k in adata.obsm.keys()]
                if getattr(adata, "isbacked", False):
                    adata.file.close()
                self.root.after(0, lambda: self._on_inspect_success(obs_columns, var_names, obsm_keys))
            except Exception as exc:
                self.root.after(0, lambda: self._on_inspect_error(exc))

        threading.Thread(target=worker, daemon=True).start()

    def _on_inspect_success(
        self, obs_columns: Sequence[str], var_names: Sequence[str], obsm_keys: Sequence[str]
    ) -> None:
        values = sorted(set(str(col) for col in obs_columns))
        self._obs_columns = list(values)
        self._var_names = sorted(set(str(g) for g in var_names))
        self._obsm_keys = sorted(set(str(k) for k in obsm_keys))

        self.section_key_combo.configure(values=self._obs_columns)
        self.color_combo.configure(values=self._obs_columns)
        self.outline_combo.configure(values=[""] + self._obs_columns)
        self.cell_annotations_editor.set_choices(self._obs_columns)
        self.neighbor_stats_annotations_editor.set_choices(self._obs_columns)
        self.genes_editor.set_choices(self._var_names)
        self.spatial_key_combo.configure(values=self._obsm_keys)

        if self.spatial_key.get() not in self._obsm_keys:
            fallback_candidates = [
                "spatial",
                "X_spatial",
                "spatial_coords",
                "X_spatial_coords",
                "Spatial",
                "spatialcoords",
            ]
            found_fallback = False
            for fb in fallback_candidates:
                if fb in self._obsm_keys:
                    self.spatial_key.set(fb)
                    found_fallback = True
                    break
            if not found_fallback and self._obsm_keys:
                self.spatial_key.set(self._obsm_keys[0])

        if self.section_key.get() not in self._obs_columns:
            if "sample_id" in self._obs_columns:
                self.section_key.set("sample_id")
            elif self._obs_columns:
                self.section_key.set(self._obs_columns[0])

        if self.annotation.get() not in self._obs_columns:
            if "leiden" in self._obs_columns:
                self.annotation.set("leiden")
            elif self._obs_columns:
                self.annotation.set(self._obs_columns[0])

        if self.outline_by.get() and self.outline_by.get() not in self._obs_columns:
            if "course" in self._obs_columns:
                self.outline_by.set("course")
            elif "condition" in self._obs_columns:
                self.outline_by.set("condition")
            else:
                self.outline_by.set("")

        self._log(f"Inspect complete: {len(self._obs_columns)} obs columns, {len(self._var_names)} genes.")
        self._set_busy(False, "Ready.")

    def _on_inspect_error(self, exc: Exception) -> None:
        assert messagebox is not None
        tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__)).strip()
        self._log("Inspect failed:")
        self._log(tb if tb else f"{type(exc).__name__}: {exc}")
        self._set_busy(False, "Inspect failed.")
        messagebox.showerror("Inspect failed", f"{type(exc).__name__}: {exc}")

    def _validate_and_collect(self) -> Tuple[str, dict, dict]:
        input_raw = self.input_path.get().strip()
        output_raw = self.output_path.get().strip()
        if not input_raw:
            raise ValueError("Input .h5ad path is required.")
        if not output_raw:
            raise ValueError("Output .html path is required.")

        input_path = Path(input_raw).expanduser().resolve()
        if not input_path.exists():
            raise ValueError(f"Input file not found: {input_path}")
        if input_path.suffix.lower() != ".h5ad":
            raise ValueError("Input file must have .h5ad extension.")

        output_path = Path(output_raw).expanduser()
        if output_path.suffix.lower() != ".html":
            output_path = output_path.with_suffix(".html")
        if not output_path.is_absolute():
            output_path = (Path.cwd() / output_path).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)

        section_key = self.section_key.get().strip()
        annotation = self.annotation.get().strip()
        if not section_key:
            raise ValueError("section_key is required.")
        if not annotation:
            raise ValueError("Main cells annotation is required.")

        section_metadata = _unique(_parse_tokens(self.section_metadata.get() or ""))
        metadata_max_columns_raw = self.metadata_max_columns.get().strip()
        metadata_max_columns = None
        if metadata_max_columns_raw:
            metadata_max_columns = _parse_non_negative_int("Metadata max columns", metadata_max_columns_raw)

        metadata_value_order_raw = self.metadata_value_order_text.get("1.0", "end").strip()
        metadata_value_order = None
        if metadata_value_order_raw and metadata_value_order_raw != "{}":
            try:
                parsed = json.loads(metadata_value_order_raw)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Metadata value order must be valid JSON: {exc}") from exc
            if not isinstance(parsed, dict):
                raise ValueError("Metadata value order JSON must be an object/dictionary.")
            normalized = {}
            for key, values in parsed.items():
                if values is None:
                    normalized[str(key)] = []
                elif isinstance(values, list):
                    normalized[str(key)] = [str(v) for v in values]
                else:
                    raise ValueError("Each metadata value order entry must be a list.")
            metadata_value_order = normalized

        metadata_labels_raw = self.metadata_labels_text.get("1.0", "end").strip()
        metadata_labels = None
        if metadata_labels_raw and metadata_labels_raw != "{}":
            try:
                parsed = json.loads(metadata_labels_raw)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Metadata labels must be valid JSON: {exc}") from exc
            if not isinstance(parsed, dict):
                raise ValueError("Metadata labels JSON must be an object/dictionary.")
            metadata_labels = {
                str(key): str(value)
                for key, value in parsed.items()
                if str(key).strip() and value is not None
            } or None

        min_panel_size = _parse_positive_int("Min panel size", self.min_panel_size.get())
        spot_size = _parse_spot_size(self.spot_size.get())

        downsample_raw = self.downsample.get().strip()
        downsample = None
        if downsample_raw:
            downsample = _parse_positive_int("Downsample", downsample_raw)

        feature_sparse_zero_threshold = float(self.feature_sparse_zero_threshold.get().strip())
        if feature_sparse_zero_threshold < 0 or feature_sparse_zero_threshold > 1:
            raise ValueError("Feature sparse threshold must be between 0 and 1.")
        feature_storage = self.feature_storage.get().strip() or "embedded"
        feature_manifest_path_raw = self.feature_manifest_path.get().strip()
        feature_manifest_path = None
        if feature_manifest_path_raw:
            feature_manifest_path = str(Path(feature_manifest_path_raw).expanduser())

        gene_correlation_top_n = _parse_non_negative_int("Corr. top N", self.gene_correlation_top_n.get())
        spatial_variable_genes_n = _parse_non_negative_int("SVG n", self.spatial_variable_genes_n.get())
        cell_annotations = _unique(self.cell_annotations_editor.get_items())
        features = _unique(self.genes_editor.get_items())
        outline_by = self.outline_by.get().strip() or None

        load_kwargs = {
            "section_key": section_key,
            "section_metadata": section_metadata,
            "metadata_value_order": metadata_value_order,
            "metadata_max_columns": metadata_max_columns,
            "spatial_key": self.spatial_key.get().strip() or "spatial",
        }

        export_kwargs = {
            "output_path": str(output_path),
            "main_cell_annotation": annotation,
            "title": self.title.get().strip() or "KaroSpace",
            "min_panel_size": min_panel_size,
            "spot_size": spot_size,
            "downsample": downsample,
            "outline_by": outline_by,
            "metadata_labels": metadata_labels,
            "cell_annotations": cell_annotations,
            "features": features,
            "gene_correlation_top_n": gene_correlation_top_n,
            "spatial_variable_genes_n": spatial_variable_genes_n,
            "feature_encoding": self.feature_encoding.get().strip() or "auto",
            "feature_storage": feature_storage,
            "feature_manifest_path": feature_manifest_path,
            "feature_sparse_zero_threshold": feature_sparse_zero_threshold,
            "neighbor_stats_permutations": _parse_neighbor_permutations(self.neighbor_permutations.get()),
            "neighbor_stats_seed": int(self.neighbor_stats_seed.get().strip() or "0"),
            "interaction_markers_top_targets": _parse_positive_int("Interaction top targets", self.interaction_markers_top_targets.get()),
            "interaction_markers_top_genes": _parse_positive_int("Interaction top genes", self.interaction_markers_top_genes.get()),
            "interaction_markers_min_cells": _parse_positive_int("Interaction min cells", self.interaction_markers_min_cells.get()),
            "interaction_markers_min_neighbors": _parse_positive_int(
                "Interaction min neighbors",
                self.interaction_markers_min_neighbors.get(),
            ),
        }

        if bool(self.neighbor_stats_auto.get()):
            export_kwargs["neighbor_stats_annotations"] = [color]
        else:
            export_kwargs["neighbor_stats_annotations"] = _unique(self.neighbor_stats_annotations_editor.get_items())

        return str(input_path), load_kwargs, export_kwargs

    def export_html(self) -> None:
        if self._is_busy:
            return
        assert messagebox is not None
        try:
            input_path, load_kwargs, export_kwargs = self._validate_and_collect()
        except Exception as exc:
            messagebox.showerror("Invalid settings", str(exc))
            return

        self._set_busy(True, "Exporting...")
        self._log(f"Loading dataset: {input_path}")
        self._log(
            f"section_key={load_kwargs['section_key']}, "
            f"main_cell_annotation={export_kwargs['main_cell_annotation']}"
        )
        if export_kwargs.get("cell_annotations"):
            self._log(f"cell_annotations={len(export_kwargs['cell_annotations'])}")
        if export_kwargs.get("features"):
            self._log(f"features={len(export_kwargs['features'])} (manual list)")
        if export_kwargs.get("feature_storage") == "sidecar":
            self._log("feature_storage=sidecar (lazy aux loading enabled)")

        def worker() -> None:
            try:
                from .data_loader import load_spatial_data
                from .exporter import export_to_html

                dataset = load_spatial_data(input_path, **load_kwargs)
                output_path = export_to_html(dataset, **export_kwargs)
                self.root.after(0, lambda: self._on_export_success(str(output_path)))
            except Exception as exc:
                self.root.after(0, lambda: self._on_export_error(exc))

        threading.Thread(target=worker, daemon=True).start()

    def _on_export_success(self, output_path: str) -> None:
        assert messagebox is not None
        self._log(f"Export finished: {output_path}")
        self._set_busy(False, "Export complete.")
        message = f"KaroSpace HTML created:\n{output_path}"
        if self.feature_storage.get().strip() == "sidecar":
            aux_path = self.feature_manifest_path.get().strip()
            if not aux_path:
                aux_path = str(Path(output_path).with_suffix(".features.json"))
            message += f"\n\nFeature sidecar:\n{aux_path}"
        messagebox.showinfo("Export complete", message)

    def _on_export_error(self, exc: Exception) -> None:
        assert messagebox is not None
        self._log("Export failed:")
        self._log(traceback.format_exc().strip())
        self._set_busy(False, "Export failed.")
        messagebox.showerror("Export failed", str(exc))


def launch_gui(initial_input: Optional[str] = None, initial_output: Optional[str] = None) -> None:
    if tk is None or ttk is None:
        raise RuntimeError(
            "Tkinter is not available in this Python environment. Install Python with Tk support "
            "or run exports via the `karospace` CLI."
        )
    root = tk.Tk()
    KaroSpaceExportGUI(root, initial_input=initial_input, initial_output=initial_output)
    root.mainloop()


def main() -> None:
    parser = argparse.ArgumentParser(description="Launch the KaroSpaceBuilder GUI.")
    parser.add_argument("--input", type=str, default=None, help="Optional .h5ad path to prefill.")
    parser.add_argument("--output", type=str, default=None, help="Optional output HTML path to prefill.")
    args = parser.parse_args()
    launch_gui(initial_input=args.input, initial_output=args.output)


if __name__ == "__main__":
    main()
