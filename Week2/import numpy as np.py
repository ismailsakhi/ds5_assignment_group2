import numpy as np
import matplotlib.pyplot as plt

def mandelbrot_smooth(xmin, xmax, ymin, ymax, width, height, max_iter=300):
    # Grid of c values
    x = np.linspace(xmin, xmax, width).reshape(1, -1)
    y = np.linspace(ymin, ymax, height).reshape(-1, 1)
    c = x + 1j*y

    z = np.zeros_like(c, dtype=np.complex128)
    escaped = np.zeros(c.shape, dtype=bool)
    smooth = np.zeros(c.shape, dtype=np.float64)

    for i in range(max_iter):
        z = z*z + c
        mag = np.abs(z)
        just_escaped = (~escaped) & (mag > 2.0)
        if np.any(just_escaped):
            smooth[just_escaped] = (i + 1) - np.log2(np.log(mag[just_escaped]))
            escaped[just_escaped] = True
        if escaped.all():
            break

    # Points that didn't escape get max_iter
    smooth[~escaped] = max_iter
    return smooth

def zoom_bounds(xmin, xmax, ymin, ymax, cx, cy, factor):
    """Return new bounds zoomed by `factor` into center (cx, cy)."""
    x_half = (xmax - xmin) / (2.0 * factor)
    y_half = (ymax - ymin) / (2.0 * factor)
    return cx - x_half, cx + x_half, cy - y_half, cy + y_half

# --- Try it ---
width = height = 800
max_iter = 300

# Full view
xmin, xmax, ymin, ymax = -2.0, 1.0, -1.5, 1.5
img = mandelbrot_smooth(xmin, xmax, ymin, ymax, width, height, max_iter)
plt.figure(figsize=(6,6))
plt.imshow(img, extent=[xmin, xmax, ymin, ymax], origin='lower', vmin=0, vmax=max_iter)
plt.title("Mandelbrot – Full View"); plt.xlabel("Re(c)"); plt.ylabel("Im(c)")
plt.colorbar(label="Smooth iterations")
plt.show()

# Zoom 1: Seahorse Valley near (-0.75, 0.10)
cx, cy = -0.75, 0.10
xmin, xmax, ymin, ymax = zoom_bounds(xmin, xmax, ymin, ymax, cx, cy, factor=8)
img = mandelbrot_smooth(xmin, xmax, ymin, ymax, width, height, max_iter)
plt.figure(figsize=(6,6))
plt.imshow(img, extent=[xmin, xmax, ymin, ymax], origin='lower', vmin=0, vmax=max_iter)
plt.title("Zoom 1 – Around (-0.75, 0.10)")
plt.xlabel("Re(c)"); plt.ylabel("Im(c)")
plt.colorbar(label="Smooth iterations")
plt.show()

# Zoom 2: deeper
xmin, xmax, ymin, ymax = zoom_bounds(xmin, xmax, ymin, ymax, cx, cy, factor=10)
img = mandelbrot_smooth(xmin, xmax, ymin, ymax, width, height, max_iter)
plt.figure(figsize=(6,6))
plt.imshow(img, extent=[xmin, xmax, ymin, ymax], origin='lower', vmin=0, vmax=max_iter)
plt.title("Zoom 2 – Deeper at (-0.75, 0.10)")
plt.xlabel("Re(c)"); plt.ylabel("Im(c)")
plt.colorbar(label="Smooth iterations")
plt.show()


width = height = 800
max_iter = 300
init_bounds = (-2.0, 1.0, -1.5, 1.5)
bounds = list(init_bounds)
zoom_factor = 2.0

fig, ax = plt.subplots(figsize=(6,6))

def render():
    ax.clear()
    img = mandelbrot_smooth(*bounds, width, height, max_iter)
    ax.imshow(img, extent=[bounds[0], bounds[1], bounds[2], bounds[3]],
              origin='lower', vmin=0, vmax=max_iter)
    ax.set_title("Click to zoom (LMB in, RMB out). Press 'r' to reset.")
    ax.set_xlabel("Re(c)"); ax.set_ylabel("Im(c)")
    plt.draw()

def onclick(event):
    if not event.inaxes:
        return
    cx, cy = event.xdata, event.ydata
    f = zoom_factor if event.button == 1 else 1/zoom_factor
    bounds[:] = zoom_bounds(*bounds, cx, cy, factor=f)
    render()

def onkey(event):
    if event.key == 'r':
        bounds[:] = init_bounds
        render()

cid = fig.canvas.mpl_connect('button_press_event', onclick)
kid = fig.canvas.mpl_connect('key_press_event', onkey)

render()
plt.show()
