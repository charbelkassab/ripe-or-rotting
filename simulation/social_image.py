"""Social media image: a real simulated odour field with fly tracks, headline and key numbers.

  python simulation/social_image.py   -> paper/social_portrait.png (1080x1350), paper/social_landscape.png (1200x627)
"""

import json

import numpy as np
from PIL import Image, ImageDraw, ImageFilter, ImageFont

from atmosphere import Plume

BG = (8, 11, 20)
INK = (238, 241, 246)
MUTED = (140, 150, 168)
RIPE = np.array([70, 150, 255], float)
ROT = np.array([255, 128, 60], float)
SURGE = (120, 235, 150)
CAST = (150, 158, 176)
SRC = [(0.0, 0.15), (0.0, -0.15)]

def font(size, weight='regular'):
    try:
        if weight == 'bold':
            return ImageFont.truetype('/System/Library/Fonts/Avenir Next.ttc', size, index=0)
        index = {'regular': 7, 'demi': 2, 'medium': 5}[weight]
        return ImageFont.truetype('/System/Library/Fonts/Avenir Next.ttc', size, index=index)
    except OSError:
        return ImageFont.truetype('DejaVuSans.ttf', size)


def odour_field(w, h, x0, x1, y0, y1, seed=21):
    p = Plume(SRC, x_max=x1 + 0.3, seed=seed)
    for _ in range(int(40 / 0.05)):
        p.step(0.05)
    xs = np.linspace(x0, x1, w)
    ys = np.linspace(y1, y0, h)
    img = np.zeros((h, w, 3))
    for k, col in enumerate([RIPE, ROT]):
        c = p.c_unit_grid(xs, ys, source=k)
        a = np.clip((np.log10(c + 1e-9) + 0.5) / 4.0, 0, 1) ** 1.3
        img += a[..., None] * col[None, None]
    return img


def render(W, H, field_box, text_layout, out):
    fx, fy, fw, fh = field_box
    X0, X1 = -0.12, 2.25
    Y0, Y1 = -(X1 - X0) * fh / fw / 2, (X1 - X0) * fh / fw / 2
    canvas = Image.new('RGB', (W, H), BG)

    low = odour_field(fw // 3, fh // 3, X0, X1, Y0, Y1)
    glow = Image.fromarray(np.clip(low, 0, 255).astype(np.uint8)).resize((fw, fh), Image.BICUBIC)
    bloom = glow.filter(ImageFilter.GaussianBlur(18))
    field = Image.blend(bloom, glow, 0.65)
    # fade the field edges into the background
    yy, xx = np.mgrid[0:fh, 0:fw]
    edge = np.minimum.reduce([xx / 60, (fw - xx) / 160, yy / 90, (fh - yy) / 90])
    alpha = Image.fromarray((np.clip(edge, 0, 1) * 255).astype(np.uint8))
    canvas.paste(field, (fx, fy), alpha)

    d = ImageDraw.Draw(canvas)

    def px(x, y):
        return fx + (x - X0) / (X1 - X0) * fw, fy + (Y1 - y) / (Y1 - Y0) * fh

    tracks = json.load(open('results/tracks.json'))
    overlay = Image.new('RGBA', (W, H), (0, 0, 0, 0))
    od = ImageDraw.Draw(overlay)
    scored = sorted(tracks, key=lambda tr: -sum(1 for r in tr if r[3] == 1))[:16]
    for tr in scored:
        a = np.array(tr)
        pts = [(px(x, y), m) for _, x, y, m in a[::2]
               if 0.14 <= x <= X1 and Y0 * 0.85 <= y <= Y1 * 0.85]
        for (p1, m1), (p2, _) in zip(pts, pts[1:]):
            if abs(p1[0] - p2[0]) + abs(p1[1] - p2[1]) > 60:
                continue
            col = SURGE if m1 == 1 else CAST
            od.line([p1, p2], fill=col + (150,), width=2)
    canvas.paste(overlay, (0, 0), overlay)
    d = ImageDraw.Draw(canvas)

    for (x, y), col, lab in zip(SRC, [RIPE, ROT], ['ripe banana', 'rotting banana']):
        cx, cy = px(x, y)
        d.ellipse([cx - 13, cy - 13, cx + 13, cy + 13], fill=tuple(int(v) for v in col), outline=INK, width=3)
        f = font(int(fw * 0.024), 'demi')
        d.text((cx + 22, cy - 16 if y > 0 else cy - 6), lab, font=f, fill=INK, stroke_width=4, stroke_fill=BG)

    # wind arrow and scale bar
    f_small = font(int(fw * 0.02), 'medium')
    ax, ay = fx + fw - 290, fy + fh - 34
    d.line([ax, ay, ax + 90, ay], fill=MUTED, width=3)
    d.polygon([(ax + 100, ay), (ax + 86, ay - 7), (ax + 86, ay + 7)], fill=MUTED)
    d.text((ax + 110, ay - 14), 'wind 0.3 m/s', font=f_small, fill=MUTED, stroke_width=3, stroke_fill=BG)
    x_a, _ = px(0.25, 0)
    x_b, _ = px(1.25, 0)
    sy = fy + fh - 34
    d.line([x_a, sy, x_b, sy], fill=MUTED, width=3)
    d.line([x_a, sy - 7, x_a, sy + 7], fill=MUTED, width=3)
    d.line([x_b, sy - 7, x_b, sy + 7], fill=MUTED, width=3)
    d.text(((x_a + x_b) / 2 - 18, sy - 40), '1 m', font=f_small, fill=MUTED, stroke_width=3, stroke_fill=BG)

    text_layout(d, canvas)
    canvas.save(out, optimize=True)
    print('wrote', out)


def portrait():
    W, H = 1080, 1350
    M = 72

    def text(d, canvas):
        d.text((M, 70), 'SIMULATION STUDY  ·  165,122-NEURON FLY CONNECTOME', font=font(22, 'demi'), fill=(255, 150, 90))
        f_h = font(70, 'bold')
        y = 112
        for line in ['Why do fruit flies', 'end up on rotting fruit?']:
            d.text((M, y), line, font=f_h, fill=INK)
            y += 84
        d.text((M, y + 18), 'I simulated one, from banana chemistry and turbulent air',
               font=font(29, 'regular'), fill=MUTED)
        d.text((M, y + 58), 'through its entire nervous system to where it lands and eats.',
               font=font(29, 'regular'), fill=MUTED)

        y0 = 1000
        stats = [('10×', 'lower smell threshold', 'for rotting banana'),
                 ('8–19 m', 'detection range for rotting', 'vs 2–6 m for ripe'),
                 ('50 / 50', 'first landings up close:', 'smell stops mattering')]
        col_w = (W - 2 * M) // 3
        for i, (big, l1, l2) in enumerate(stats):
            x = M + i * col_w
            if i:
                d.line([x - 18, y0 + 6, x - 18, y0 + 150], fill=(40, 48, 66), width=2)
            d.text((x, y0), big, font=font(62, 'bold'), fill=INK)
            d.text((x, y0 + 84), l1, font=font(24, 'demi'), fill=(205, 212, 224))
            d.text((x, y0 + 116), l2, font=font(24, 'regular'), fill=MUTED)

        d.line([M, 1205, W - M, 1205], fill=(40, 48, 66), width=2)
        d.text((M, 1232), 'Rotting fruit is louder. The yeast on it tastes better. Hunger decides.',
               font=font(28, 'demi'), fill=(255, 150, 90))
        d.text((M, 1280), 'Charbel Kassab  ·  charbelk.com/ripe-or-rotting', font=font(24, 'medium'), fill=MUTED)
        lx, ly = M, 950
        for lab, col in [('simulated fly: searching', CAST), ('flying upwind after its brain detects the smell', SURGE)]:
            d.line([lx, ly + 13, lx + 26, ly + 13], fill=col, width=3)
            d.text((lx + 36, ly), lab, font=font(21, 'medium'), fill=MUTED)
            lx += 70 + d.textlength(lab, font=font(21, 'medium'))

    render(W, H, (0, 390, W, 550), text, 'paper/social_portrait.png')


def landscape():
    W, H = 1200, 627
    M = 56

    def text(d, canvas):
        d.text((M, 46), 'SIMULATION STUDY  ·  165,122-NEURON FLY CONNECTOME', font=font(17, 'demi'), fill=(255, 150, 90))
        f_h = font(46, 'bold')
        d.text((M, 76), 'Why do fruit flies end up on rotting fruit?', font=f_h, fill=INK)
        d.text((M, 140), 'Rotting fruit is louder. The yeast on it tastes better. Hunger decides.',
               font=font(24, 'demi'), fill=(205, 212, 224))
        d.text((M, H - 44), 'Charbel Kassab  ·  charbelk.com/ripe-or-rotting', font=font(19, 'medium'), fill=MUTED)

    render(W, H, (0, 185, W, 380), text, 'paper/social_landscape.png')


if __name__ == '__main__':
    portrait()
    landscape()
