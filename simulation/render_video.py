"""Video of the stimulus-reaction simulation: two bananas, turbulent plumes, flies surging and casting.

  .venv/bin/python render_video.py [state]   -> results/foraging_<state>.mp4 and .gif
"""

import sys

import imageio.v2 as imageio
import numpy as np
from PIL import Image, ImageDraw, ImageFont

from navigate import simulate

W, H = 1280, 720
X0, X1, Y0, Y1 = -0.15, 2.25, -0.68, 0.68
FONT = '/System/Library/Fonts/HelveticaNeue.ttc'  # macOS; falls back to DejaVu / PIL default elsewhere
BG = np.array([12, 15, 22], float)
RIPE = np.array([80, 160, 255], float)
ROT = np.array([255, 140, 60], float)
MODE_COL = {0: (170, 176, 190), 1: (120, 230, 140), 2: (255, 230, 90), 3: (90, 90, 100)}


def load_font(size, bold=False):
    for path, index in [(FONT, 1 if bold else 0), ('DejaVuSans-Bold.ttf' if bold else 'DejaVuSans.ttf', 0)]:
        try:
            return ImageFont.truetype(path, size, index=index)
        except OSError:
            continue
    return ImageFont.load_default(size)


def to_px(x, y):
    return (x - X0) / (X1 - X0) * W, (Y1 - y) / (Y1 - Y0) * H


def main():
    state = sys.argv[1] if len(sys.argv) > 1 else 'starved_24h'
    foods, positions = ['ripe', 'rotting'], [(0.0, 0.15), (0.0, -0.15)]
    res, _, frames = simulate(foods, positions, state=state, n_flies=60, seed=4, record_frames=True, t_max=45)
    print(res)
    xs = np.linspace(X0, X1, W // 4)
    ys = np.linspace(Y1, Y0, H // 4)
    GX, GY = np.meshgrid(xs, ys)
    pts = np.c_[GX.ravel(), GY.ravel()]
    f_big, f_small = load_font(30, True), load_font(20)
    out = f'results/foraging_{state}.mp4'
    writer = imageio.get_writer(out, fps=20, quality=8, macro_block_size=1)
    gif_frames = []
    for fi, fr in enumerate(frames):
        img = np.tile(BG, (len(ys), len(xs), 1))
        if len(fr['puffs']):
            s = fr['sig']
            m = 1.0 / 25.0
            d2 = ((pts[:, None, :] - fr['puffs'][None, :, :]) ** 2).sum(-1)
            contrib = m / ((2 * np.pi) ** 1.5 * s ** 3)[None, :] * np.exp(-d2 / (2 * s * s)[None, :])
            for k, col in enumerate([RIPE, ROT]):
                c = contrib[:, fr['sid'] == k].sum(1).reshape(GX.shape)
                a = np.clip(np.log10(c + 1e-9) / 4.0, 0, 1)[..., None] ** 1.5
                img = img * (1 - 0.85 * a) + col * 0.85 * a
        im = Image.fromarray(img.astype(np.uint8)).resize((W, H), Image.BICUBIC)
        d = ImageDraw.Draw(im)
        for (x, y), col, lab in zip(positions, [RIPE, ROT], ['ripe banana', 'rotting banana (day 7)']):
            px, py = to_px(x, y)
            d.ellipse([px - 14, py - 14, px + 14, py + 14], fill=tuple(int(v) for v in col), outline=(255, 255, 255), width=2)
            d.text((px + 22, py - 12), lab, font=f_small, fill=(235, 235, 240))
        for (x, y), m in zip(fr['flies'], fr['mode']):
            if m == 3:
                continue
            px, py = to_px(x, y)
            r = 5 if m < 2 else 7
            d.ellipse([px - r, py - r, px + r, py + r], fill=MODE_COL[int(m)])
        fed = fr['fed_on']
        d.text((24, 18), f'Flies searching for banana in a 0.3 m/s breeze  ·  state: {state.replace("_", " ")}',
               font=f_big, fill=(240, 240, 245))
        d.text((24, 62), f't = {fr["t"]:4.1f} s   fed on ripe: {(fed == 0).sum()}   fed on rotting: {(fed == 1).sum()}',
               font=f_small, fill=(200, 205, 215))
        lx = 24
        for lab, col in [('casting (no odour)', MODE_COL[0]), ('surging upwind (odour detected by the brain)', MODE_COL[1]),
                         ('feeding', MODE_COL[2])]:
            d.ellipse([lx, H - 38, lx + 14, H - 24], fill=col)
            d.text((lx + 22, H - 42), lab, font=f_small, fill=(200, 205, 215))
            lx += 48 + d.textlength(lab, font=f_small)
        d.text((W - 330, H - 42), 'wind  →   (2× speed)', font=f_small, fill=(200, 205, 215))
        arr = np.asarray(im)
        writer.append_data(arr)  # one frame per 0.1 s of simulation at 20 fps -> 2x speed
        if fi % 3 == 0:
            gif_frames.append(im.resize((640, 360)))
    writer.close()
    gif_frames[0].save(f'results/foraging_{state}.gif', save_all=True, append_images=gif_frames[1:], duration=150, loop=0)
    print('wrote', out)


if __name__ == '__main__':
    main()
