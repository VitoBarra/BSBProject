import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def tukey_summary(values: np.ndarray) -> dict:
    """Return the quartiles and whisker limits used by a standard box plot."""
    q1, median, q3 = np.percentile(values, [25, 50, 75])
    iqr = q3 - q1
    lower_fence = q1 - 1.5 * iqr
    upper_fence = q3 + 1.5 * iqr

    inlier_values = values[(values >= lower_fence) & (values <= upper_fence)]
    whisker_low = np.min(inlier_values)
    whisker_high = np.max(inlier_values)

    return {
        "q1": q1,
        "median": median,
        "q3": q3,
        "iqr": iqr,
        "lower_fence": lower_fence,
        "upper_fence": upper_fence,
        "whisker_low": whisker_low,
        "whisker_high": whisker_high,
    }


def plot_distribution(ax_data, ax_box, values: np.ndarray, title: str, color: str) -> None:
    """Show the ordered data with quantile guides and the related box plot."""
    sorted_values = np.sort(values)
    x_positions = np.arange(1, len(sorted_values) + 1)
    stats = tukey_summary(sorted_values)

    ax_data.plot(
        x_positions,
        sorted_values,
        marker="o",
        linestyle="-",
        linewidth=1,
        markersize=4,
        color=color,
        alpha=0.9,
    )

    guide_lines = [
        ("Q1", stats["q1"], "#1f77b4"),
        ("Median", stats["median"], "#d62728"),
        ("Q3", stats["q3"], "#2ca02c"),
        ("Whisker low", stats["whisker_low"], "#7f7f7f"),
        ("Whisker high", stats["whisker_high"], "#7f7f7f"),
    ]

    for label, y_value, line_color in guide_lines:
        ax_data.axhline(
            y_value,
            color=line_color,
            linestyle="--",
            linewidth=1,
            label=f"{label}: {y_value:.1f}",
        )

    ax_data.set_title(f"{title}\nOrdered discrete data")
    ax_data.set_ylabel("Value")
    ax_data.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
    ax_data.set_xlim(1, len(sorted_values))
    ax_data.legend(loc="upper left", fontsize=8, frameon=True)

    ax_box.boxplot(
        sorted_values,
        vert=True,
        patch_artist=True,
        boxprops={"facecolor": color, "alpha": 0.5},
        medianprops={"color": "#d62728", "linewidth": 2},
        whiskerprops={"color": "#444444"},
        capprops={"color": "#444444"},
        flierprops={
            "marker": "o",
            "markerfacecolor": "#111111",
            "markeredgecolor": "#111111",
            "markersize": 5,
            "alpha": 0.8,
        },
    )
    ax_box.set_title("Matching box plot")
    ax_box.set_ylabel("Value")
    ax_box.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.5)


rng = np.random.default_rng(7)

# Each dataset is intentionally small enough that the ordered observations stay readable.
datasets = [
    (
        "Nearly symmetric",
        np.round(rng.normal(loc=50, scale=7, size=24)).astype(int),
        "#4c78a8",
    ),
    (
        "Right-skewed",
        np.round(rng.lognormal(mean=2.2, sigma=0.45, size=24)).astype(int),
        "#f58518",
    ),
    (
        "Bimodal",
        np.sort(
            np.concatenate(
                [
                    np.round(rng.normal(loc=25, scale=3, size=12)),
                    np.round(rng.normal(loc=55, scale=4, size=12)),
                ]
            )
        ).astype(int),
        "#54a24b",
    ),
    (
        "Compact with outliers",
        np.array([18, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 22, 23, 23, 24, 24, 25, 26, 8, 34]),
        "#e45756",
    ),
]

fig, axes = plt.subplots(
    nrows=len(datasets),
    ncols=2,
    figsize=(12, 16),
    gridspec_kw={"width_ratios": [4, 1.2]},
)

for row_index, (name, values, color) in enumerate(datasets):
    plot_distribution(axes[row_index, 0], axes[row_index, 1], values, name, color)

fig.suptitle(
    "How ordered data, quartiles, whiskers, and box plots relate",
    fontsize=16,
    y=0.995,
)
fig.supxlabel("Ordered observation", fontsize=12)
plt.tight_layout(rect=(0, 0, 1, 0.985))

output_path = Path(__file__).with_name("boxplot_gallery.png")
fig.savefig(output_path, dpi=200)
print(f"Saved figure to: {output_path}")

plt.show()
