"""Generate printable field-kit PDFs: 14-day count sheet and A4 scale bar."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

OUT = Path(__file__).resolve().parent
plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Droid Sans Fallback"]
plt.rcParams["axes.unicode_minus"] = False


def count_sheet():
    fig = plt.figure(figsize=(210 / 25.4, 297 / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 210); ax.set_ylim(0, 297); ax.axis("off")

    ax.text(105, 285, "楼梯人流计数表 / Stair Traffic Count Sheet",
            ha="center", va="center", fontsize=15, weight="bold")
    ax.text(20, 272, "建筑 / Building：____________________", fontsize=9)
    ax.text(115, 272, "楼梯编号 / Stair ID：__________", fontsize=9)
    ax.text(20, 263, "计数人 / Counter：__________________", fontsize=9)
    ax.text(115, 263, "日期范围 / Period：__________", fontsize=9)

    headers = ["日期\nDate", "星期\nDay", "08:00–08:30\n上 Up", "08:00–08:30\n下 Down",
               "12:00–12:30\n上 Up", "12:00–12:30\n下 Down", "17:30–18:00\n上 Up",
               "17:30–18:00\n下 Down", "备注\nNotes"]
    x0, y0 = 15, 245
    col_w = [22, 14, 19, 19, 19, 19, 19, 19, 30]
    row_h = 12.5
    # header
    x = x0
    for h, w in zip(headers, col_w):
        ax.add_patch(Rectangle((x, y0), w, row_h, fill=False, lw=0.8))
        ax.text(x + w / 2, y0 + row_h / 2, h, ha="center", va="center", fontsize=7)
        x += w
    # 14 rows
    for i in range(14):
        y = y0 - (i + 1) * row_h
        x = x0
        for j, w in enumerate(col_w):
            ax.add_patch(Rectangle((x, y), w, row_h, fill=False, lw=0.5))
            if j == 0:
                ax.text(x + 3, y + row_h / 2, "____ / ____", va="center", fontsize=6.5)
            if j == 1:
                ax.text(x + w / 2, y + row_h / 2, "____", ha="center", va="center", fontsize=6.5)
            x += w
    ax.text(15, 38,
            "填写说明 / Instructions\n"
            "1. 每个时段连续计数 30 分钟；同一天三个时段尽量固定。\n"
            "2. 上行/下行分开计数；无法区分时只填总数并在备注说明。\n"
            "3. 闭馆、施工、暴雨等异常日请在备注注明；不要补填估算值。\n"
            "4. 连续 14 天；若某天缺失，顺延补足 14 天。\n"
            "5. 计数时不记录任何个人身份信息。",
            fontsize=8, va="top", linespacing=1.6)
    ax.text(15, 12, "打印：A4，100%（不要缩放）", fontsize=7, color="0.4")
    fig.savefig(OUT / "count_sheet_A4.pdf")
    fig.savefig(OUT / "count_sheet_A4.png", dpi=200)
    plt.close(fig)
    print("wrote count_sheet_A4.pdf")


def scale_bar():
    fig = plt.figure(figsize=(210 / 25.4, 297 / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 210); ax.set_ylim(0, 297); ax.axis("off")

    ax.text(105, 285, "摄影测量比例尺 / Photogrammetry Scale Bar",
            ha="center", va="center", fontsize=15, weight="bold")
    ax.text(105, 276, "打印时必须选择 100%（实际大小），打印后用直尺核对 100 mm 刻度",
            ha="center", va="center", fontsize=8.5, color="0.25")

    # 100 mm bar in 10 mm blocks, alternating black/white, 10 mm tall
    x0, y0 = 55, 240
    for i in range(10):
        color = "black" if i % 2 == 0 else "white"
        ax.add_patch(Rectangle((x0 + i * 10, y0), 10, 10, facecolor=color,
                               edgecolor="black", lw=0.6))
    for i in range(11):
        ax.plot([x0 + i * 10, x0 + i * 10], [y0 - 2, y0 + 12], color="black", lw=0.6)
        if i % 5 == 0:
            ax.text(x0 + i * 10, y0 - 5, f"{i * 10}", ha="center", va="top", fontsize=7)
    ax.text(x0 + 50, y0 + 15, "100 mm", ha="center", fontsize=9)

    # 20 mm check square
    ax.add_patch(Rectangle((20, 240), 20, 20, facecolor="none", edgecolor="black", lw=1))
    ax.text(30, 235, "20 mm 核对", ha="center", fontsize=7)

    # 50 x 50 mm checkerboard (5 mm squares)
    bx, by = 80, 150
    for i in range(10):
        for j in range(10):
            color = "black" if (i + j) % 2 == 0 else "white"
            ax.add_patch(Rectangle((bx + i * 5, by + j * 5), 5, 5,
                                   facecolor=color, edgecolor="black", lw=0.2))
    ax.text(bx + 25, by - 6, "50 × 50 mm 棋盘格 / checkerboard", ha="center", fontsize=7)

    # coded targets: simple black squares with white centres
    for k, (cx, cy) in enumerate([(30, 170), (150, 170), (30, 200), (150, 200)]):
        ax.add_patch(Rectangle((cx, cy), 14, 14, facecolor="black", edgecolor="black"))
        ax.add_patch(Rectangle((cx + 4, cy + 4), 6, 6, facecolor="white", edgecolor="white"))
        ax.text(cx + 7, cy - 3, f"T{k+1}", ha="center", fontsize=6)

    ax.text(105, 120,
            "使用说明 / How to use\n"
            "1. 把本页平放在踏面上，确保不变形、不反光。\n"
            "2. 至少两个比例尺，互相垂直；棋盘格用于检查镜头畸变。\n"
            "3. 拍摄时让比例尺出现在每一组照片中；不要裁掉。\n"
            "4. 每组 40–80 张，重叠 70%；避免正午直射光。\n"
            "5. 打印后先量 100 mm 刻度；若误差 >1 mm，重新打印。",
            fontsize=8.5, va="top", linespacing=1.7)
    fig.savefig(OUT / "scale_bar_A4.pdf")
    fig.savefig(OUT / "scale_bar_A4.png", dpi=200)
    plt.close(fig)
    print("wrote scale_bar_A4.pdf")


if __name__ == "__main__":
    count_sheet()
    scale_bar()
