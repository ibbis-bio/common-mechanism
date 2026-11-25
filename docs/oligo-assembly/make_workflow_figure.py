"""Annotate the existing workflow figure with a plain-language findings legend.

Keeps the original flowchart on the left. Places a readable legend on the right
with a vertical separator. Does not push to GitHub.
"""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent
BACKUP = ROOT / "workflow-original-backup.png"
OUT = ROOT / "workflow.png"


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
    return lines


def bullet_block(draw, x, y, title, items, title_font, body_font, max_width, gap=8):
    for line in wrap(draw, title, title_font, max_width):
        draw.text((x, y), line, fill=(35, 40, 90), font=title_font)
        y += 24
    y += 6
    for item in items:
        lines = wrap(draw, item, body_font, max_width - 20)
        draw.text((x, y), "•", fill=(35, 40, 90), font=body_font)
        for i, line in enumerate(lines):
            draw.text((x + 16, y + i * 19), line, fill=(60, 65, 86), font=body_font)
        y += max(1, len(lines)) * 19 + gap
    return y


def main():
    if not BACKUP.exists():
        raise SystemExit(f"Missing original backup at {BACKUP}")

    base = Image.open(BACKUP).convert("RGBA")
    left_w, height = base.size

    title_f = font(22, bold=True)
    head_f = font(16, bold=True)
    body_f = font(13)
    intro_f = font(13)

    panel_w = 620
    gap = 28
    sep_x = left_w + gap // 2
    panel_x = left_w + gap
    total_w = left_w + gap + panel_w

    intro = (
        "An oligo is a short DNA fragment in the pool. COMMEC is the screening tool that looks "
        "for sequences of concern. OliGraph does not screen. It only tries to join short fragments "
        "into longer assembled pieces called contigs. COMMEC then screens those contigs. The report "
        "shows one row for each original oligo."
    )
    mains = [
        "COMMEC direct: COMMEC found a concern in this short fragment itself, before any assembly.",
        "COMMEC on OliGraph contig: this fragment alone did not show a concern, but COMMEC found a concern after OliGraph joined it into a longer piece.",
        "COMMEC direct and OliGraph contig: COMMEC found a concern both in the original fragment and in the longer assembled piece.",
        "No concern: COMMEC did not find a concern in the original fragment or in any assembled piece that includes it.",
        "Incomplete: screening could not finish for a practical reason, such as the fragment being too short, too long, not assemblable, or a failed screen.",
    ]
    details = [
        "2A: COMMEC found a concern in this oligo. Step 2A means the direct screen of the original fragment.",
        "4B: COMMEC found a concern in an OliGraph contig. Step 2B built the contig. Step 4B screened it.",
        "2A+4B concern: COMMEC found a concern in both the original fragment and its assembled contig.",
        "2A+4B no concern: COMMEC checked both the original fragment and its contig and found no concern.",
        "2A skipped because the oligo was too short. 4B found no concern in the assembled contig.",
        "2A and 4B both skipped because the oligo and the contig were too short to screen.",
        "2B: OliGraph built no contig that includes this oligo, so there was no assembled piece to screen.",
        "2B: this oligo was not eligible for OliGraph assembly.",
        "4B: COMMEC started screening the assembled contig but did not finish.",
        "2A: COMMEC failed while screening this original oligo.",
    ]

    # Measure content height so the legend can be vertically centered.
    probe = Image.new("RGBA", (10, 10), (255, 255, 255, 255))
    probe_draw = ImageDraw.Draw(probe)
    content_width = panel_w - 40
    content_h = 0
    content_h += 32  # title
    content_h += len(wrap(probe_draw, intro, intro_f, content_width)) * 18
    content_h += 14
    content_h += 30  # main heading
    for item in mains:
        content_h += max(1, len(wrap(probe_draw, item, body_f, content_width - 20))) * 19 + 8
    content_h += 10 + 30  # detail heading
    for item in details:
        content_h += max(1, len(wrap(probe_draw, item, body_f, content_width - 20))) * 19 + 8

    total_h = max(height, content_h + 48)
    canvas = Image.new("RGBA", (total_w, total_h), (255, 255, 255, 255))
    # Center the flowchart vertically if the legend is taller.
    top = max(0, (total_h - height) // 2)
    canvas.paste(base, (0, top))
    draw = ImageDraw.Draw(canvas)

    # Vertical separator between flowchart and legend.
    draw.line([(sep_x, 20), (sep_x, total_h - 20)], fill=(35, 40, 90), width=3)
    draw.line([(sep_x + 4, 20), (sep_x + 4, total_h - 20)], fill=(180, 186, 200), width=1)

    # Legend panel border.
    panel_top = 16
    panel_bottom = total_h - 16
    draw.rectangle(
        [panel_x, panel_top, total_w - 16, panel_bottom],
        fill=(248, 249, 252, 255),
        outline=(35, 40, 90, 255),
        width=2,
    )

    x = panel_x + 18
    max_w = content_width
    panel_inner_h = panel_bottom - panel_top
    y = panel_top + max(12, (panel_inner_h - content_h) // 2)

    draw.text((x, y), "How to read the findings", fill=(35, 40, 90), font=title_f)
    y += 32
    for line in wrap(draw, intro, intro_f, max_w):
        draw.text((x, y), line, fill=(60, 65, 86), font=intro_f)
        y += 18

    y += 14
    y = bullet_block(
        draw,
        x,
        y,
        "Main findings summary",
        mains,
        head_f,
        body_f,
        max_w,
        gap=8,
    )
    y += 10
    bullet_block(
        draw,
        x,
        y,
        "Detailed findings",
        details,
        head_f,
        body_f,
        max_w,
        gap=8,
    )

    canvas.convert("RGB").save(OUT, format="PNG", optimize=True)
    print(f"Wrote side-by-side annotated workflow to {OUT} ({total_w}x{total_h})")


if __name__ == "__main__":
    main()
