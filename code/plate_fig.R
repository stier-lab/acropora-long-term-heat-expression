# =============================================================================
# Purpose: RNA-seq 96-well plate-layout figure (Figure B) — a batch-robust
#          plating design balancing treatment, wound, genet, and day across
#          plates, with extraction/PCR controls and anchor samples marked.
#          Molly's original layout script, ported into the repo for provenance.
# Input:   none (layout is defined inline)
# Output:  figures/Fig_B_balanced_plates.{pdf,png}  (if `fig_dir` is set,
#          else written to the working directory)
#          figures/Fig_B_plate_layouts.csv          — the well-by-well layout
# Note:    STANDALONE — not sourced by code/_run_all.R (run directly when the
#          sequencing plate map is needed). Uses its own palette, independent
#          of 00_setup.R, to match the wet-lab plate-map conventions.
# =============================================================================

# Requires: tidyverse, patchwork
library(tidyverse)
library(patchwork)

# -----------------------------
# 1) Palettes & helpers
# -----------------------------
pal <- list(
  temp   = c(`28`="#1A9BD7", `31`="#FF6F5B"),  # cool vs warm
  ctrl   = "#B3A739",                          # controls (olive)
  anchor = "#36B39C",                          # anchors (teal)
  spare  = "#E6A8F7",                          # spare (pink)
  wound  = c(U="#6B6B6B", W="#151515"),        # border: U gray; W black
  other  = "#8A8A8A",                          # border for non-primaries
  spare_border = "#C9C9C9"                     # border for spare wells
)

make_plate_grid <- function(plate_id){
  expand.grid(row = LETTERS[1:8], col = 1:12) %>%
    mutate(plate = plate_id, x = col, y = 9 - pmin(match(row, LETTERS), 8))
}

# -----------------------------
# 2) Controls & anchor positions
# -----------------------------
controls_pos <- tribble(
  ~row, ~col, ~control_type,
  "A", 1,  "Extraction NTC",
  "A", 2,  "Extraction NTC",
  "B", 1,  "Library NTC",
  "B", 12, "Library NTC",
  "H", 11, "Spike-in",
  "H", 12, "Spike-in",
  "G", 1,  "Tech Rep (anchor)",
  "G", 12, "Tech Rep (anchor)"
)

# No rename needed — tribble already names columns row/col
anchors_pos <- tribble(
  ~row, ~col,
  "D", 6,
  "D", 7,
  "E", 6,
  "E", 7,
  "C", 5,
  "C", 8,
  "F", 5,
  "F", 8
)

# Safety checks (bounds & overlaps)
.valid_rows <- LETTERS[1:8]
stopifnot(all(controls_pos$row %in% .valid_rows), all(anchors_pos$row %in% .valid_rows))
stopifnot(all(controls_pos$col %in% 1:12),       all(anchors_pos$col %in% 1:12))
.overlap <- inner_join(controls_pos %>% select(row,col),
                       anchors_pos  %>% select(row,col), by = c("row","col"))
if(nrow(.overlap)>0) stop("Control and anchor wells overlap at: ",
                          paste0(.overlap$row,.overlap$col, collapse=", "))

# -----------------------------
# 3) Experimental design (balanced)
# -----------------------------
temps <- c(28,31); days <- c(0,3,9); wounds <- c("U","W")
tanks <- c("A","B","C","D"); genotypes <- c("A","B","C")

design <- expand.grid(Temp=temps, Day=days, Wound=wounds, Tank=tanks, Geno=genotypes,
                      stringsAsFactors = FALSE) %>%
  arrange(Temp, Day, Wound, Tank, Geno)

set.seed(123)
design <- design %>%
  group_by(Temp, Day, Wound) %>%
  mutate(Plate = sample(rep(c(1,2), length.out = n()))) %>%
  ungroup()

# -----------------------------
# 4) Reserve wells & fill plates
# -----------------------------
mark_reserved <- function(plate_df, include_anchors=TRUE){
  pd <- plate_df %>%
    left_join(controls_pos, by=c("row","col")) %>%
    mutate(is_control = !is.na(control_type))
  if(include_anchors){
    pd <- pd %>%
      left_join(anchors_pos %>% mutate(is_anchor=TRUE), by=c("row","col")) %>%
      mutate(is_anchor = replace_na(is_anchor, FALSE))
  } else pd <- pd %>% mutate(is_anchor=FALSE)
  pd
}

fill_plate <- function(plate_df, sample_tbl){
  avail <- plate_df %>% filter(!is_control, !is_anchor) %>% arrange(row, col)
  n_assign <- min(nrow(sample_tbl), nrow(avail))
  primaries <- bind_cols(avail[seq_len(n_assign),], sample_tbl[seq_len(n_assign),]) %>%
    mutate(SampleType="Primary", control_type=NA_character_)
  controls <- plate_df %>% filter(is_control) %>%
    mutate(SampleType="Control", Temp=NA_real_, Day=NA_real_,
           Wound=NA_character_, Tank=NA_character_, Geno=NA_character_,
           control_type = control_type %||% NA_character_)
  anchors <- plate_df %>% filter(is_anchor) %>%
    mutate(SampleType="Anchor", Temp=NA_real_, Day=NA_real_,
           Wound=NA_character_, Tank=NA_character_, Geno=NA_character_,
           control_type="Anchor")
  used <- bind_rows(controls %>% select(row,col),
                    anchors %>% select(row,col),
                    primaries %>% select(row,col))
  spares <- anti_join(plate_df, used, by=c("row","col")) %>%
    mutate(SampleType="Spare", Temp=NA_real_, Day=NA_real_,
           Wound=NA_character_, Tank=NA_character_, Geno=NA_character_,
           control_type=NA_character_)
  common_cols <- c("row","col","plate","x","y","is_control","control_type","is_anchor",
                   "SampleType","Temp","Day","Wound","Tank","Geno")
  normalize <- function(df){
    df %>% mutate(
      Temp=as.numeric(Temp), Day=as.numeric(Day),
      Wound=as.character(Wound), Tank=as.character(Tank), Geno=as.character(Geno),
      control_type=as.character(control_type), SampleType=as.character(SampleType)
    ) %>% select(all_of(common_cols)) %>% mutate(well=paste0(row,col))
  }
  bind_rows(normalize(controls), normalize(anchors),
            normalize(primaries), normalize(spares))
}

plate1 <- mark_reserved(make_plate_grid(1), include_anchors=TRUE)
plate2 <- mark_reserved(make_plate_grid(2), include_anchors=TRUE)
samples_p1 <- design %>% filter(Plate==1) %>% arrange(Temp, Day, Wound, Tank, Geno)
samples_p2 <- design %>% filter(Plate==2) %>% arrange(Temp, Day, Wound, Tank, Geno)

plate1_filled <- fill_plate(plate1, samples_p1)
plate2_filled <- fill_plate(plate2, samples_p2)

# -----------------------------
# 5) Build plotting data
# -----------------------------
plates_all <- bind_rows(plate1_filled, plate2_filled) %>%
  mutate(
    Plate   = factor(plate, levels=c(1,2), labels=c("Plate 1","Plate 2")),
    TankNum = rep(1:12, length.out=n()),
    Thicket = rep(c("A","C","D"), length.out=n()),
    lbl_top = if_else(SampleType=="Primary", paste0("D", Day, " ", Wound), ""),
    lbl_bot = case_when(
      SampleType != "Primary" ~ "",
      TRUE ~ paste0("T", TankNum, "·", Thicket)
    ),
    ctrl_lbl = recode(control_type,
                      "Extraction NTC"     = "Extr NTC",
                      "Library NTC"        = "Lib NTC",
                      "Spike-in"           = "Spike-in",
                      "Tech Rep (anchor)"  = "Tech Rep",
                      .default = control_type),
    stroke_col = case_when(
      SampleType=="Primary" & Wound=="U" ~ pal$wound["U"],
      SampleType=="Primary" & Wound=="W" ~ pal$wound["W"],
      SampleType=="Control"              ~ pal$other,
      SampleType=="Anchor"               ~ pal$other,
      TRUE                               ~ pal$spare_border
    ),
    stroke_lwd = case_when(
      SampleType=="Primary" ~ 0.8,
      SampleType=="Control" ~ 0.6,
      SampleType=="Anchor"  ~ 0.6,
      TRUE                  ~ 0.5
    ),
    fill_val = case_when(
      SampleType=="Primary" ~ paste0("T", Temp),
      SampleType=="Control" ~ "Control",
      SampleType=="Anchor"  ~ "Anchor",
      TRUE                  ~ "Spare"
    ),
    lbl_top_sh = lbl_top,
    lbl_bot_sh = lbl_bot
  )

# Fail-fast balance check
bad <- plates_all %>%
  filter(SampleType=="Primary") %>%
  count(Temp, Day, Wound, plate, name="n") %>%
  pivot_wider(names_from=plate, values_from=n, values_fill=0) %>%
  filter(`1`!=6 | `2`!=6)
if(nrow(bad)>0) stop("Imbalance in Temp×Day×Wound:\n", paste(capture.output(print(bad)), collapse="\n"))

# -----------------------------
# 6) Plot
# -----------------------------
base <- ggplot(plates_all, aes(x=x,y=y)) +
  geom_tile(aes(fill = fill_val), width=.95, height=.95, color=NA) +
  geom_rect(aes(xmin=x-.475, xmax=x+.475, ymin=y-.475, ymax=y+.475),
            color = plates_all$stroke_col, linewidth = plates_all$stroke_lwd, fill=NA) +
  geom_text(aes(label=lbl_top_sh), vjust=0.0, nudge_y= 0.15, size=2.9, color=NA) +
  geom_text(aes(label=lbl_bot_sh), vjust=1.2, nudge_y=-0.15, size=2.9, color=NA) +
  geom_text(aes(label=lbl_top),    vjust=0.0, nudge_y= 0.15, size=2.8) +
  geom_text(aes(label=lbl_bot),    vjust=1.2, nudge_y=-0.15, size=2.8) +
  geom_text(
    data = subset(plates_all, SampleType %in% c("Control","Anchor")),
    aes(label = ifelse(SampleType=="Control", ctrl_lbl, "Anchor")),
    size = 2.6, lineheight = 0.95
  ) +
  scale_x_continuous(breaks=1:12, expand=expansion(add=.3)) +
  scale_y_continuous(breaks=1:8 , expand=expansion(add=.3)) +
  facet_wrap(~ Plate, nrow=1) + coord_equal() +
  scale_fill_manual(
    name=NULL,
    breaks=c("T28","T31","Control","Anchor","Spare"),
    labels=c("28 °C (Cool)","31 °C (Warm)","Control","Anchor","Spare"),
    values=c(T28=unname(pal$temp["28"]), T31=unname(pal$temp["31"]),
             Control=pal$ctrl, Anchor=pal$anchor, Spare=pal$spare)
  ) +
  labs(
    title="Balanced Plating & Controls (2 × 96-well plates)",
    subtitle="Fill: Temperature (28 °C blue; 31 °C coral). Border: Wound (U gray; W black). Controls/Anchors labeled.",
    x="Column", y="Row"
  ) +
  theme_minimal(base_size=11) +
  theme(
    panel.grid=element_blank(), axis.ticks=element_blank(),
    strip.text=element_text(face="bold", size=12),
    plot.title=element_text(face="bold", size=15),
    plot.subtitle=element_text(size=10, color="#4a4a4a"),
    legend.position="none"
  )

# Combined legend row: fill + wound-border
legend_fill <- tibble(
  Category = factor(c("28 °C (Cool)", "31 °C (Warm)", "Control", "Anchor", "Spare"),
                    levels = c("28 °C (Cool)", "31 °C (Warm)", "Control", "Anchor", "Spare")),
  Color = c(unname(pal$temp["28"]), unname(pal$temp["31"]), pal$ctrl, pal$anchor, pal$spare),
  x = seq(1, 15, by = 3),  # increase horizontal spacing
  y = 1
) %>%
  ggplot(aes(x, y, fill = Color)) +
  geom_tile(width = 2.5, height = 0.9, color = "black") +  # wider tiles
  geom_text(aes(label = Category), size = 4, fontface = "bold", vjust = 0.5, hjust = 0.5) +
  scale_fill_identity() +
  coord_equal(xlim = c(0, 17)) +
  theme_void() +
  theme(
    plot.margin = margin(0, 8, 0, 0),
    plot.background = element_rect(fill = "white", color = NA)
  )

legend_wound <- tibble(Wound = c("U", "W"), x = c(1, 4), y = 1) %>%
  ggplot(aes(x, y)) +
  # Draw wider tiles with borders
  geom_tile(width = 2.5, height = 1.0, fill = "#FFFFFF", color = "black") +
  geom_rect(
    aes(xmin = x - 1.25, xmax = x + 1.25, ymin = y - 0.5, ymax = y + 0.5),
    color = c(pal$wound["U"], pal$wound["W"]),
    linewidth = 1.4,
    fill = NA
  ) +
  # Place text INSIDE each cell, centered
  annotate("text", x = 1, y = 1, label = "Unwounded border", size = 3.8, fontface = "bold", vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 4, y = 1, label = "Wounded border", size = 3.8, fontface = "bold", vjust = 0.5, hjust = 0.5) +
  coord_equal(xlim = c(0, 6)) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 8),
    plot.background = element_rect(fill = "white", color = NA)
  )


final_fig <- base / (legend_fill | legend_wound) + plot_layout(heights=c(10,1))
print(final_fig)

# Save
out_pdf <- if (exists("fig_dir")) file.path(fig_dir, "Fig_B_balanced_plates.pdf") else "Fig_B_balanced_plates.pdf"
out_png <- if (exists("fig_dir")) file.path(fig_dir, "Fig_B_balanced_plates.png") else "Fig_B_balanced_plates.png"

# High-res PDF (vector) and PNG (600 dpi)
ggsave(out_pdf, final_fig, width = 12, height = 7, dpi = 600, bg = "white")

# Prefer ragg device for the PNG if available (crisper text rendering)
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(
    filename   = out_png,
    width      = 20,
    height     = 7,
    units      = "in",
    res        = 600,
    background = "white"
  )
  print(final_fig)
  grDevices::dev.off()
} else {
  ggsave(out_png, final_fig, width = 12, height = 7, dpi = 600, bg = "white")
}

# Optional: export plate layouts as CSVs
plate_layouts <- plates_all %>%
  filter(SampleType=="Primary") %>%
  transmute(Plate, well, Temp, Day, Wound, Tank, Geno)
write_csv(plate_layouts, if (exists("fig_dir")) file.path(fig_dir,"Fig_B_plate_layouts.csv") else "Fig_B_plate_layouts.csv")