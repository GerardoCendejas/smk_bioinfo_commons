#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

option_list <- list(
  make_option(c("--ind1"), type="character", help="TSV for Population 1"),
  make_option(c("--ind2"), type="character", help="TSV for Population 2"),
  make_option(c("--ind3"), type="character", help="TSV for Population 3"),
  make_option(c("--comb12"), type="character", help="TSV for Combined 1+2"),
  make_option(c("--comb13"), type="character", help="TSV for Combined 1+3"),
  make_option(c("--comb23"), type="character", help="TSV for Combined 2+3"),
  make_option(c("--sp1"), type="character", default="Pop1"),
  make_option(c("--sp2"), type="character", default="Pop2"),
  make_option(c("--sp3"), type="character", default="Pop3"),
  make_option(c("-o", "--output"), type="character", default="ccr_plot.png"),
  make_option(c("-g", "--generation_time"), type="numeric", default=1.0),
  make_option(c("-s", "--start_generations"), type = "numeric", default = 0, 
              help = "Filter out generations before this value [default = %default]", metavar = "number"),
  make_option(c("-c", "--confidence"), type = "numeric", default = 0.95,
              help = "Confidence interval for CCR [default = %default]", metavar = "number"),
  make_option(c("--gif_output"), type="character", default=NULL,
              help="Output GIF path for Ne/CCR animation (requires gganimate + gifski) [default: none]"),
  make_option(c("--gif3d_output"), type="character", default=NULL,
              help="Output GIF path for 3D rotating population tube plot (requires gganimate + gifski) [default: none]")
)

args <- parse_args(OptionParser(option_list=option_list))

load_ne <- function(path) {
    df <- read.delim(path, header=TRUE, check.names=FALSE)
    time <- df[[1]]
    ne_matrix <- df[, -1, drop=FALSE]
    return(list(time=time, ne=as.matrix(ne_matrix)))
}

load_eta <- function(path) {
    df <- read.delim(path, header=TRUE, check.names=FALSE)
    time <- df[[1]]
    eta_matrix <- 1 / (2 * df[, -1, drop=FALSE])
    return(list(time=time, eta=as.matrix(eta_matrix)))
}

# Each Phlash run is independent and has its own time slices, will use step_interp to evaluate at
# a common grid T, this is in a similar manner as in the phlash paper @ https://github.com/jthlab/phlash_paper/blob/master/notebooks/ccr.py
step_interp <- function(time, values, T) {
  approx(time, values, xout = T, method = "constant", rule = 2)$y
}

eta_on_grid <- function(time, eta_matrix, T) {
  apply(eta_matrix, 2, function(col) step_interp(time, col, T))
}

calculate_ccr <- function(time_a, eta_a, time_b, eta_b, time_comb, eta_comb, T, label, gen_time) {
  eta_a_grid    <- eta_on_grid(time_a, eta_a, T)
  eta_b_grid    <- eta_on_grid(time_b, eta_b, T)
  eta_comb_grid <- eta_on_grid(time_comb, eta_comb, T)
  ccr_matrix <- (2 * eta_comb_grid) / (eta_a_grid + eta_b_grid)

  data.frame(
    Time = T * gen_time,
    Median = apply(ccr_matrix, 1, median),
    Low = apply(ccr_matrix, 1, quantile, probs=(1 - args$confidence) / 2),
    High = apply(ccr_matrix, 1, quantile, probs=1 - (1 - args$confidence) / 2),
    Comparison = label
  )
}

pop1 <- load_eta(args$ind1)
pop2 <- load_eta(args$ind2)
pop3 <- load_eta(args$ind3)
comb12 <- load_eta(args$comb12)
comb13 <- load_eta(args$comb13)
comb23 <- load_eta(args$comb23)

# Shared grid, so the times are actually present in all the files, no extrapolation outside any actual inference
all_times <- list(pop1$time, pop2$time, pop3$time, comb12$time, comb13$time, comb23$time)
t_min <- max(sapply(all_times, min))
t_max <- min(sapply(all_times, max))
T_grid <- exp(seq(log(t_min), log(t_max), length.out = 1000))

pair12 <- calculate_ccr(pop1$time, pop1$eta, pop2$time, pop2$eta, comb12$time, comb12$eta, T_grid,
                        paste(args$sp1, "vs", args$sp2), args$generation_time)

pair13 <- calculate_ccr(pop1$time, pop1$eta, pop3$time, pop3$eta, comb13$time, comb13$eta, T_grid,
                        paste(args$sp1, "vs", args$sp3), args$generation_time)

pair23 <- calculate_ccr(pop2$time, pop2$eta, pop3$time, pop3$eta, comb23$time, comb23$eta, T_grid,
                        paste(args$sp2, "vs", args$sp3), args$generation_time)

plot_data <- rbind(pair12, pair13, pair23)

plot_data <- plot_data %>% filter(Time >= args$start_generations * args$generation_time)

p <- ggplot(plot_data, aes(x = Time, y = Median, linetype = Comparison)) +
  geom_ribbon(aes(ymin = Low, ymax = High), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) + 
  labs(
    x = ifelse(args$generation_time == 1, "Time (Generations)", "Time (Years)"),
    y = "Relative Cross-Coalescent Rate",
    title = "Population Divergence History",
    subtitle = expression(frac(eta[combined], eta[1] + eta[2]))
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(face="bold"))

ggsave(args$output, plot = p, width = 10, height = 7, dpi = 300)
cat("Success! CCR plot saved to:", args$output, "\n")

needs_anim <- !is.null(args$gif_output) || !is.null(args$gif3d_output)
if (needs_anim) {
  if (!requireNamespace("gganimate", quietly=TRUE) || !requireNamespace("gifski", quietly=TRUE))
    stop("--gif_output / --gif3d_output requires R packages: gganimate, gifski")
  suppressPackageStartupMessages(library(gganimate))

  # CCR → distance: d = (1−CCR)/CCR; CCR=1→d=0 (panmictic), CCR=0.5→d=1 (threshold), CCR→0→d=∞ (fully diverged)
  # ponytail: -log(CCR) matches PNG log-scale axis and avoids 1/CCR explosion for small CCR; CCR>1 clamped to 0
  ccr_df <- data.frame(
    Time = pair12$Time,
    d12  = pmax(0, -log(pair12$Median)),
    d13  = pmax(0, -log(pair13$Median)),
    d23  = pmax(0, -log(pair23$Median))
  )

  pop1_ne <- load_ne(args$ind1)
  pop2_ne <- load_ne(args$ind2)
  pop3_ne <- load_ne(args$ind3)

# Also use the common grid
  ne_df <- data.frame(
    Time = T_grid * args$generation_time,
    ne1  = apply(eta_on_grid(pop1_ne$time, pop1_ne$ne, T_grid), 1, median),
    ne2  = apply(eta_on_grid(pop2_ne$time, pop2_ne$ne, T_grid), 1, median),
    ne3  = apply(eta_on_grid(pop3_ne$time, pop3_ne$ne, T_grid), 1, median)
  )

  combined <- merge(ccr_df, ne_df, by = "Time")
  combined  <- combined[combined$Time >= args$start_generations * args$generation_time, ]

  # Classical MDS via cmdscale (base R): exact 2D embedding when triangle inequality holds,
  # degrades gracefully (discards negative eigenvalues) when CCR distances violate it —
  # which trilaterate could not handle (it degenerated to a collinear configuration with
  # one point flung far along the x-axis, producing the spurious extreme distances).
  embed_3pts <- function(d12, d13, d23) {
    D   <- matrix(c(0, d12, d13, d12, 0, d23, d13, d23, 0), nrow = 3)
    pos <- suppressWarnings(cmdscale(D, k = 2))
    # cmdscale returns fewer than 2 cols when <2 positive eigenvalues (coincident or collinear points)
    if (!is.matrix(pos) || ncol(pos) < 2) {
      pad <- matrix(0, 3, 2)
      if (is.matrix(pos) && ncol(pos) == 1) pad[, 1] <- pos[, 1]
      return(pad)
    }
    pos
  }

  # Sort in animation order (ancient→recent) for sequential Procrustes alignment
  combined <- combined[order(combined$Time, decreasing = TRUE), ]

  anim_rows <- lapply(seq_len(nrow(combined)), function(i) {
    pos <- embed_3pts(combined$d12[i], combined$d13[i], combined$d23[i])
    data.frame(Time = combined$Time[i],
               pop  = c(args$sp1, args$sp2, args$sp3),
               x    = pos[, 1], y = pos[, 2],
               Ne   = c(combined$ne1[i], combined$ne2[i], combined$ne3[i]))
  })

  # Sequential no-reflection Procrustes: prevents left-right orientation flips between frames.
  # Each frame is rotated (never reflected) to best match the previous frame's configuration.
  positions <- lapply(anim_rows, function(df) as.matrix(df[, c("x", "y")]))
  for (i in seq(2, length(positions))) {
    sv  <- svd(t(positions[[i]]) %*% positions[[i - 1]])   # P_new^T P_ref
    d   <- sign(det(sv$u %*% t(sv$v)))                     # +1 = pure rotation, no mirror
    R   <- sv$u %*% diag(c(1, d)) %*% t(sv$v)
    positions[[i]] <- positions[[i]] %*% R
  }

  anim_data <- do.call(rbind, lapply(seq_along(anim_rows), function(i) {
    df <- anim_rows[[i]]; df$x <- positions[[i]][, 1]; df$y <- positions[[i]][, 2]; df
  }))
  # ponytail: log10(Time) so each decade gets equal animation time; distant past flies, recent crawls
  anim_data$t_anim <- -log10(anim_data$Time)
  time_unit        <- ifelse(args$generation_time == 1, "generations", "years")

  if (!is.null(args$gif_output)) {
    # Fixed global limits so distances reflect absolute divergence magnitude across frames
    xpad  <- diff(range(anim_data$x)) * 0.1
    ypad  <- diff(range(anim_data$y)) * 0.1
    xlims <- range(anim_data$x) + c(-xpad, xpad)
    ylims <- range(anim_data$y) + c(-ypad, ypad)

    p_anim <- ggplot(anim_data, aes(x = x, y = y, color = pop)) +
      geom_point(aes(size = Ne), alpha = 0.85) +
      scale_color_brewer(palette = "Set1") +
      scale_size_continuous(trans = "log10", range = c(3, 20)) +
      scale_x_continuous(limits = xlims) +
      scale_y_continuous(limits = ylims) +
      labs(
        title  = paste0("{formatC(10^(-frame_time), format='e', digits=2)} ", time_unit, " ago"),
        x      = NULL, y = NULL,
        color  = "Population", size = "Ne"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom", legend.box = "horizontal") +
      transition_time(t_anim) +
      ease_aes("linear")

    animate(p_anim, nframes = 150, fps = 12, width = 700, height = 650,
            renderer = gifski_renderer(args$gif_output))
    cat("Animation saved to:", args$gif_output, "\n")
  }

  if (!is.null(args$gif3d_output)) {
    # 3D tubes: x,y = trilaterate coords, z = log10(time) normalized to [0,1] (0=present, 1=past)
    traj     <- anim_data[order(anim_data$pop, anim_data$Time), ]
    z_raw    <- log10(traj$Time)
    traj$z   <- (z_raw - min(z_raw)) / diff(range(z_raw))
    # Normalize x,y with same scale to preserve pairwise CCR distances in the projection
    xy_scale <- max(diff(range(traj$x)), diff(range(traj$y)), 1e-9)
    traj$xn  <- (traj$x - mean(range(traj$x))) / xy_scale
    traj$yn  <- (traj$y - mean(range(traj$y))) / xy_scale

    phi    <- pi / 4   # elevation: 45° so z-axis and xy-plane share equal visual weight
    n_rot  <- 90
    thetas <- seq(0, 2 * pi, length.out = n_rot + 1)[seq_len(n_rot)]

    # Orthographic turntable projection: rotate around z, then tilt by phi
    rot_df <- do.call(rbind, lapply(seq_along(thetas), function(fi) {
      theta <- thetas[fi]
      xr    <- traj$xn * cos(theta) - traj$yn * sin(theta)
      yr    <- traj$xn * sin(theta) + traj$yn * cos(theta)
      data.frame(
        Time  = traj$Time, pop = traj$pop, Ne = traj$Ne,
        xs    = xr,
        ys    = -yr * sin(phi) + traj$z * cos(phi),
        frame = fi
      )
    }))

    p3d <- ggplot(rot_df, aes(x = xs, y = ys, color = pop, group = pop)) +
      geom_path(linewidth = 1.5) +
      geom_point(aes(size = Ne), alpha = 0.7) +
      scale_color_brewer(palette = "Set1") +
      scale_size_continuous(trans = "log10", range = c(2, 10)) +
      coord_fixed() +
      labs(
        title  = paste0("Population divergence trajectories  |  bottom = present, top = past (log ", time_unit, ")"),
        color  = "Population", size = "Ne", x = NULL, y = NULL
      ) +
      theme_void(base_size = 14) +
      theme(legend.position = "bottom", legend.box = "horizontal",
            plot.title = element_text(face = "bold", hjust = 0.5, size = 11)) +
      transition_manual(frame)

    animate(p3d, nframes = n_rot, fps = 15, width = 700, height = 700,
            renderer = gifski_renderer(args$gif3d_output))
    cat("3D rotation GIF saved to:", args$gif3d_output, "\n")
  }
}
