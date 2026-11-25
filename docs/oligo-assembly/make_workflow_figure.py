"""Annotate the existing workflow figure with HTML-matched findings badges.

Keeps the original flowchart on the left. Places a legend on the right that uses
the same Main Finding colors and short badge labels as the HTML report.
Does not push to GitHub.
"""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent
BACKUP = ROOT / "workflow-original-backup.png"
OUT = ROOT / "workflow.png"

# Match report.js findingMainColor / findingMainBadge.
MAIN_COLORS = {
    "COMMEC direct": (240, 80, 35),
    "COMMEC on OliGraph contig": (107, 76, 154),
    "COMMEC direct and OliGraph contig": (194, 58, 20),
    "No concern": (65, 155, 185),
    "Incomplete": (138, 143, 158),
}
MAIN_BADGES = {
    "COMMEC direct": "Direct",
    "COMMEC on OliGraph contig": "Contig",
    "COMMEC direct and OliGraph contig": "Direct + contig",
    "No concern": "No concern",
    "Incomplete": "Incomplete",
}

MAIN_ITEMS = [
    (
        "COMMEC direct",
        "COMMEC found a concern in this short fragment itself, before any assembly. Workflow step 2A.",
    ),
    (
        "COMMEC on OliGraph contig",
        "This fragment alone did not show a concern. COMMEC found a concern after OliGraph joined it into a longer contig. Steps 2B then 4B.",
    ),
    (
        "COMMEC direct and OliGraph contig",
        "COMMEC found a concern both in the original fragment and in the longer assembled piece. Steps 2A and 4B.",
    ),
    (
        "No concern",
        "COMMEC did not find a concern in the original fragment or in any assembled contig that includes it. Steps 2A and 4B.",
    ),
    (
        "Incomplete",
        "Screening could not finish for a practical reason, such as too short, too long, no contig, exclusion, or a failed screen.",
    ),
]

# Detailed findings grouped like the HTML How to read trees.
DETAIL_GROUPS = [
    (
        "2A",
        [
            ("COMMEC direct", "COMMEC found a concern in this oligo. Step 2A is the direct screen of the original fragment."),
            ("Incomplete", "Skipped because the oligo was too short. 4B found no concern in the assembled contig."),
            ("Incomplete", "COMMEC failed while screening this original oligo."),
        ],
    ),
    (
        "2B",
        [
            ("Incomplete", "OliGraph built no contig that includes this oligo, so there was no assembled piece to screen."),
            ("Incomplete", "This oligo was not eligible for OliGraph assembly."),
        ],
    ),
    (
        "4B",
        [
            ("COMMEC on OliGraph contig", "COMMEC found a concern in an OliGraph contig. Step 2B built the contig. Step 4B screened it."),
            ("Incomplete", "COMMEC started screening the assembled contig but did not finish."),
        ],
    ),
    (
        "2A+4B",
        [
            ("COMMEC direct and OliGraph contig", "COMMEC found a concern in both the original fragment and its assembled contig."),
            ("No concern", "COMMEC checked both the original fragment and its contig and found no concern."),
            ("Incomplete", "Both skipped because the oligo and the contig were too short to screen."),
        ],
    ),
]


def font(size, bold=False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for path in candidates:
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            continue
    return ImageFont.load_default()


def wrap(draw, text, font_obj, max_width):
    words = text.split()
    lines = []
    current = ""
    for word in words:
        trial = word if not current else current + " " + word
        if draw.textlength(trial, font=font_obj) <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    return lines or [""]


def rounded_rect(draw, box, radius, fill, outline=None, width=1):
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def draw_badge(draw, x, y, label, color, badge_font, width=132, height=28):
    rounded_rect(draw, [x, y, x + width, y + height], 999, fill=color)
    tw = draw.textlength(label, font=badge_font)
    tx = x + max(4, (width - tw) / 2)
    ty = y + (height - 14) / 2
    draw.text((tx, ty), label, fill=(255, 255, 255), font=badge_font)
    return height


def draw_step_pill(draw, x, y, label, pill_font, width=72, height=28):
    rounded_rect(draw, [x, y, x + width, y + height], 6, fill=(35, 40, 90))
    tw = draw.textlength(label, font=pill_font)
    draw.text((x + (width - tw) / 2, y + (height - 13) / 2), label, fill=(255, 255, 255), font=pill_font)
    return height


def main():
    if not BACKUP.exists():
        raise SystemExit(f"Missing original backup at {BACKUP}")

    base = Image.open(BACKUP).convert("RGBA")
    left_w, height = base.size

    title_f = font(22, bold=True)
    head_f = font(15, bold=True)
    body_f = font(12)
    intro_f = font(12)
    badge_f = font(11, bold=True)
    pill_f = font(11, bold=True)
    mini_f = font(10, bold=True)

    panel_w = 640
    gap = 28
    sep_x = left_w + gap // 2
    panel_x = left_w + gap
    total_w = left_w + gap + panel_w
    content_width = panel_w - 44
    badge_w = 132

    intro = (
        "An oligo is a short DNA fragment in the pool. Status is the final combined call after step 5. "
        "Main Finding says whether a concern came from COMMEC on the original fragment or from COMMEC on "
        "an OliGraph contig. Detailed findings name the workflow steps. OliGraph only joins fragments. "
        "COMMEC always does the screening."
    )

    # Probe content height.
    probe = Image.new("RGBA", (10, 10), (255, 255, 255, 255))
    pd = ImageDraw.Draw(probe)
    content_h = 28  # title
    content_h += len(wrap(pd, intro, intro_f, content_width)) * 16 + 12
    content_h += 26  # Main heading
    for _, text in MAIN_ITEMS:
        lines = wrap(pd, text, body_f, content_width - badge_w - 16)
        content_h += max(32, len(lines) * 16) + 10
    content_h += 18 + 26  # Detail heading
    for _, items in DETAIL_GROUPS:
        content_h += 8
        for _, text in items:
            lines = wrap(pd, text, body_f, content_width - 72 - 28 - 16)
            content_h += max(40, 18 + len(lines) * 15) + 8

    total_h = max(height, content_h + 56)
    canvas = Image.new("RGBA", (total_w, total_h), (255, 255, 255, 255))
    top = max(0, (total_h - height) // 2)
    canvas.paste(base, (0, top))
    draw = ImageDraw.Draw(canvas)

    draw.line([(sep_x, 20), (sep_x, total_h - 20)], fill=(35, 40, 90), width=3)
    draw.line([(sep_x + 4, 20), (sep_x + 4, total_h - 20)], fill=(180, 186, 200), width=1)

    panel_top = 16
    panel_bottom = total_h - 16
    rounded_rect(
        draw,
        [panel_x, panel_top, total_w - 16, panel_bottom],
        10,
        fill=(248, 249, 252, 255),
        outline=(35, 40, 90, 255),
        width=2,
    )

    x = panel_x + 20
    y = panel_top + 18
    draw.text((x, y), "How to read the findings", fill=(35, 40, 90), font=title_f)
    y += 28
    for line in wrap(draw, intro, intro_f, content_width):
        draw.text((x, y), line, fill=(76, 82, 98), font=intro_f)
        y += 16
    y += 12

    # Cap style section title with underline.
    draw.text((x, y), "MAIN FINDINGS SUMMARY", fill=(35, 40, 90), font=head_f)
    y += 20
    draw.line([(x, y), (x + content_width, y)], fill=(228, 231, 238), width=1)
    y += 12

    for main, text in MAIN_ITEMS:
        color = MAIN_COLORS[main]
        label = MAIN_BADGES[main]
        text_x = x + badge_w + 14
        lines = wrap(draw, text, body_f, content_width - badge_w - 16)
        block_h = max(32, len(lines) * 16)
        draw_badge(draw, x, y + max(0, (block_h - 28) // 2), label, color, badge_f, width=badge_w)
        ty = y + max(0, (block_h - len(lines) * 16) // 2)
        for line in lines:
            draw.text((text_x, ty), line, fill=(60, 65, 86), font=body_f)
            ty += 16
        y += block_h + 10

    y += 8
    draw.text((x, y), "DETAILED FINDINGS", fill=(35, 40, 90), font=head_f)
    y += 20
    draw.line([(x, y), (x + content_width, y)], fill=(228, 231, 238), width=1)
    y += 12

    for step, items in DETAIL_GROUPS:
        # Step label column + branch cards.
        group_top = y
        step_x = x
        card_x = x + 80
        card_w = content_width - 80
        for idx, (main, text) in enumerate(items):
            color = MAIN_COLORS[main]
            short = MAIN_BADGES[main]
            lines = wrap(draw, text, body_f, card_w - 24)
            card_h = 18 + max(1, len(lines)) * 15 + 12
            # Tree rail.
            rail_x = card_x - 18
            is_last = idx == len(items) - 1
            if idx == 0:
                draw_step_pill(draw, step_x, y + 4, step, pill_f, width=68, height=26)
            # vertical rail
            rail_top = y + 8
            rail_bottom = y + card_h - (card_h // 2 if is_last else 0)
            draw.line([(rail_x, rail_top), (rail_x, rail_bottom)], fill=color + (180,), width=2)
            draw.line([(rail_x, y + 14), (rail_x + 12, y + 14)], fill=color + (180,), width=2)
            draw.ellipse([rail_x - 5, y + 9, rail_x + 5, y + 19], fill=color, outline=(255, 255, 255), width=2)

            # Card with colored left border.
            rounded_rect(
                draw,
                [card_x, y, card_x + card_w, y + card_h],
                8,
                fill=(255, 255, 255, 255),
                outline=(228, 231, 238, 255),
                width=1,
            )
            draw.rectangle([card_x, y + 2, card_x + 4, y + card_h - 2], fill=color)
            # Mini badge pill
            bw = max(54, int(draw.textlength(short, font=mini_f) + 16))
            rounded_rect(draw, [card_x + 12, y + 8, card_x + 12 + bw, y + 24], 999, fill=color)
            draw.text((card_x + 20, y + 10), short, fill=(255, 255, 255), font=mini_f)
            ty = y + 28
            for line in lines:
                draw.text((card_x + 12, ty), line, fill=(60, 65, 86), font=body_f)
                ty += 15
            y += card_h + 8
        y = max(y, group_top + 36) + 6

    canvas.convert("RGB").save(OUT, format="PNG", optimize=True)
    print(f"Wrote HTML-matched findings legend to {OUT} ({total_w}x{total_h})")


if __name__ == "__main__":
    main()
