"""Build DEVELOPMENT_PROTOCOL.pdf with all formulas rendered in LaTeX.

Strategy
--------
- reportlab Platypus for the document skeleton.
- Every mathematical expression (inline AND display) is rendered with
  matplotlib's mathtext (LaTeX syntax) into a transparent PNG and
  embedded.
- No Unicode subscripts/superscripts appear anywhere -> no black-box
  glyphs.
"""

import io
import os
import hashlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image as PILImage

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT, TA_CENTER
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle,
    Preformatted,
)

OUT_PDF = r"D:\AI\MagneticField_OneCylinder\DEVELOPMENT_PROTOCOL.pdf"
import os as _os
if _os.path.exists(OUT_PDF):
    try:
        _os.remove(OUT_PDF)
    except PermissionError:
        OUT_PDF = OUT_PDF.replace(".pdf", "_v2.pdf")
TMP_DIR = r"D:\AI\MagneticField_OneCylinder\_pdfmath"
os.makedirs(TMP_DIR, exist_ok=True)

DPI = 220


# ------------------------------------------------------------------ math
def _render(expr, fontsize, dpi=DPI):
    key = hashlib.md5(f"{expr}|{fontsize}|{dpi}".encode("utf-8")).hexdigest()
    path = os.path.join(TMP_DIR, f"m_{key}.png")
    if not os.path.exists(path):
        fig = plt.figure(figsize=(8, 1))
        fig.patch.set_alpha(0.0)
        fig.text(0.01, 0.5, f"${expr}$", fontsize=fontsize,
                 va="center", ha="left", color="black")
        fig.savefig(path, format="png", dpi=dpi,
                    bbox_inches="tight", pad_inches=0.02, transparent=True)
        plt.close(fig)
    img = PILImage.open(path)
    w_pt = img.size[0] / dpi * 72.0
    h_pt = img.size[1] / dpi * 72.0
    return path, w_pt, h_pt


def im(expr, fontsize=10):
    """Inline LaTeX math -- returns an HTML <img/> snippet for use
    inside Paragraph text."""
    path, w, h = _render(expr, fontsize=fontsize)
    # Slight downscale so the image sits comfortably on a 14-pt line.
    s = 0.85
    w *= s; h *= s
    p = path.replace("\\", "/")
    return f'<img src="{p}" width="{w:.2f}" height="{h:.2f}" valign="-2"/>'


def block(expr, fontsize=14, max_width=None):
    """Display LaTeX math -- returns a centered Image flowable."""
    path, w, h = _render(expr, fontsize=fontsize)
    if max_width and w > max_width:
        s = max_width / w
        w *= s; h *= s
    img = Image(path, width=w, height=h)
    img.hAlign = "CENTER"
    return img


# ------------------------------------------------------------------ styles
styles   = getSampleStyleSheet()
H1       = ParagraphStyle("H1", parent=styles["Heading1"],
                          spaceBefore=14, spaceAfter=8,
                          textColor=colors.HexColor("#1a3a6c"))
H2       = ParagraphStyle("H2", parent=styles["Heading2"],
                          spaceBefore=10, spaceAfter=6,
                          textColor=colors.HexColor("#1a3a6c"))
BODY     = ParagraphStyle("BODY", parent=styles["BodyText"],
                          fontSize=10.5, leading=15, alignment=TA_LEFT,
                          spaceAfter=6)
CELL     = ParagraphStyle("CELL", parent=BODY,
                          fontSize=9.5, leading=13, spaceAfter=0)
CODE     = ParagraphStyle("CODE", parent=styles["Code"],
                          fontSize=9, leading=11, leftIndent=12,
                          textColor=colors.HexColor("#2b2b2b"),
                          backColor=colors.HexColor("#f4f4f4"),
                          borderPadding=4)
TITLE    = ParagraphStyle("TITLE", parent=styles["Title"],
                          fontSize=20, leading=24, alignment=TA_CENTER,
                          textColor=colors.HexColor("#1a3a6c"),
                          spaceAfter=4)
SUBTITLE = ParagraphStyle("SUBTITLE", parent=styles["Title"],
                          fontSize=13, leading=16, alignment=TA_CENTER,
                          textColor=colors.HexColor("#444444"),
                          spaceAfter=18)


def C(*args):
    """Concatenate Paragraph fragments / strings."""
    return Paragraph("".join(str(a) for a in args), BODY)


def Cell(*args):
    return Paragraph("".join(str(a) for a in args), CELL)


# ------------------------------------------------------------------ build
doc = SimpleDocTemplate(
    OUT_PDF, pagesize=A4,
    leftMargin=2.0*cm, rightMargin=2.0*cm,
    topMargin=1.8*cm, bottomMargin=1.8*cm,
    title="Development Protocol -- Cylinder Magnetic Field",
)
story = []

# --- title
story += [
    Paragraph("Development Protocol", TITLE),
    Paragraph("Magnetic field and gradient of an axially magnetized cylinder",
              SUBTITLE),
    C("This document records the engineering protocol followed to develop "
      "the two MATLAB programs in this directory:"),
    C("&bull; <b>cylinder_magnet.m</b> &mdash; brute-force surface integral "
      "via <font face='Courier'>integral2</font>"),
    C("&bull; <b>cylinder_magnet_analytic.m</b> &mdash; analytical "
      , im(r"\phi'"), "-integration plus 1-D radial Gauss&ndash;Legendre "
      "quadrature"),
]

# --- 1
story += [
    Paragraph("1. Problem statement", H1),
    C("Compute the magnetic field ", im(r"\mathbf{B}(x,y,z)"),
      " and its full gradient tensor ",
      im(r"\partial B_\alpha/\partial x_\beta"),
      " outside a uniformly axially magnetized cylinder of radius ",
      im(r"R"), ", length ", im(r"L"), " and magnetization ",
      im(r"\mathbf{M}=M_0\,\hat{z}"), " with ", im(r"\mu_0 M_0 = 1\,\mathrm{T}"),
      ", axis aligned with the Z-axis."),
    C("Required outputs:"),
    C("&bull; 1-D on-axis profile of ", im(r"B_z"), " for validation;"),
    C("&bull; 2-D maps of ", im(r"B_x"), ", ", im(r"B_y"), ", ", im(r"B_z"),
      " and the gradient tensor in an XY plane at a chosen height ",
      im(r"z = z_0"), " above the magnet;"),
    C("&bull; visualization of in-plane field direction and gradient direction."),
]

# --- 2
story += [
    Paragraph("2. Physics formulation", H1),
    Paragraph("2.1 Scalar magnetic potential", H2),
    C("For uniform magnetization the body is replaced by surface magnetic "
      "charges with density ",
      im(r"\sigma_m=\mathbf{M}\!\cdot\!\hat{\mathbf{n}}"),
      ". Only the two flat end faces carry charge for axial magnetization "
      "(top: ", im(r"z'=+L/2,\;\sigma=+M_0"),
      "; bottom: ", im(r"z'=-L/2,\;\sigma=-M_0"), ")."),
    block(r"\varphi_m(\mathbf{r})\;=\;\frac{1}{4\pi}\,"
          r"\iint_{S}\frac{\sigma_m(\mathbf{r}')}{|\mathbf{r}-\mathbf{r}'|}\,dS'"),
    Spacer(1, 4),
    Paragraph("2.2 Field and gradient", H2),
    C("Outside the magnet ", im(r"\mathbf{H}=-\nabla\varphi_m"), " and ",
      im(r"\mathbf{B}=\mu_0\mathbf{H}"), ". Both ", im(r"\mathbf{B}"), " and ",
      im(r"\partial B_\alpha/\partial x_\beta"),
      " are obtained by differentiating <i>under the integral sign</i> "
      "analytically:"),
    block(r"H_\alpha\;=\;\frac{\sigma}{4\pi}\iint_{S}"
          r"\frac{\alpha-\alpha'}{d^{3}}\,dS'"),
    block(r"\frac{\partial H_\alpha}{\partial x_\beta}\;=\;"
          r"\frac{\sigma}{4\pi}\iint_{S}\!\left[\,"
          r"-\frac{\delta_{\alpha\beta}}{d^{3}}+"
          r"\frac{3(\alpha-\alpha')(\beta-\beta')}{d^{5}}\,\right]dS'"),
    Spacer(1, 4),
    Paragraph("2.3 Validation references", H2),
    C("<b>Near / intermediate field</b>: closed form on the cylinder axis "
      "(", im(r"\rho=0"), "):"),
    block(r"B_z(0,0,z)\;=\;\frac{\mu_0 M_0}{2}\!\left["
          r"\frac{z+L/2}{\sqrt{R^{2}+(z+L/2)^{2}}}-"
          r"\frac{z-L/2}{\sqrt{R^{2}+(z-L/2)^{2}}}\right]"),
    Spacer(1, 4),
    C("<b>Far field</b>: point-dipole formula with moment ",
      im(r"\mathbf{m}=M_0\,\pi R^{2} L\,\hat{z}"), ","),
    block(r"\mathbf{B}_{\mathrm{dip}}(\mathbf{r})\;=\;"
          r"\frac{\mu_0}{4\pi}\,"
          r"\frac{3(\mathbf{m}\!\cdot\!\hat{r})\hat{r}-\mathbf{m}}{r^{3}}"),
]

# --- 3
story += [
    Paragraph("3. Implementation strategy", H1),
    C("Two implementations were produced, each a logical step of the same plan."),
    Paragraph("3.1 Brute-force (cylinder_magnet.m)", H2),
    C("&bull; Each of the three field components and each of the five "
      "independent gradient-tensor components is computed by "
      "<font face='Courier'>integral2</font> per observation point on the "
      "(x, y) grid."),
    C("&bull; Total: ", im(r"8\times N^{2}"), " adaptive 2-D integrations."),
    C("&bull; Pros: minimal physics derivation, robust."),
    C("&bull; Cons: minutes for a ", im(r"50\times 50"),
      " grid; not useful for design iteration."),

    Paragraph("3.2 Analytical reduction (cylinder_magnet_analytic.m)", H2),
    C("<b>Step 1 &mdash; analytical ", im(r"\phi'"), " integration.</b> "
      "The ", im(r"2\pi"), " integral over the disc parameter ", im(r"\phi'"),
      " is evaluated in closed form; the result is expressed in terms of "
      "the complete elliptic integrals ", im(r"K(m)"), " and ", im(r"E(m)"),
      " with parameter"),
    block(r"m(r',\rho,\zeta)\;=\;\frac{4\,r'\,\rho}{(r'+\rho)^{2}+\zeta^{2}}"),
    Spacer(1, 4),
    C("<b>Step 2 &mdash; closed-form integrands.</b> The remaining "
      "integrand depends only on ", im(r"r',\rho,\zeta"), ":"),
    block(r"H_z(\rho,\zeta)\;=\;\frac{\sigma\,\zeta}{\pi}\,"
          r"\int_{0}^{R}\frac{r'\,E(m)}"
          r"{[(r'-\rho)^{2}+\zeta^{2}]\,\sqrt{(r'+\rho)^{2}+\zeta^{2}}}\,dr'"),
    block(r"H_\rho(\rho,\zeta)\;=\;\frac{\sigma}{2\pi\rho}\,"
          r"\int_{0}^{R}\!\!\frac{r'\,[(\rho^{2}-r'^{2}-\zeta^{2})\,E(m)+"
          r"((r'-\rho)^{2}+\zeta^{2})\,K(m)]}"
          r"{[(r'-\rho)^{2}+\zeta^{2}]\,\sqrt{(r'+\rho)^{2}+\zeta^{2}}}\,dr'"),
    Spacer(1, 4),
    C("<b>Step 3 &mdash; vectorized 1-D quadrature.</b> The radial integral "
      "on ", im(r"[0,R]"), " is evaluated by a single 96-point "
      "Gauss&ndash;Legendre rule (Golub&ndash;Welsch), broadcast across the "
      "entire observation grid as one matrix product."),
    C("<b>Step 4 &mdash; gradient tensor.</b> Six perturbation evaluations "
      "of the analytic field at observation points ", im(r"x\pm h"), ", ",
      im(r"y\pm h"), ", ", im(r"z\pm h"), " give the gradient tensor by "
      "central differences. The field is quadrature-noiseless, so a step ",
      im(r"h=10^{-6}\,\mathrm{m}"), " gives ", im(r"\sim 10^{-9}"),
      " accuracy."),
    C("Speed-up vs brute force: ", im(r"\approx 10^{3}"), "."),
]

# --- 4 : numerical-precision table (with LaTeX images in cells)
story += [Paragraph("4. Numerical-precision protocol", H1)]
issues = [
    [Cell("<b>Issue</b>"), Cell("<b>Mitigation</b>")],
    [Cell(im(r"K(m)"), " logarithmic singularity at ", im(r"m\to 1")),
     Cell("clamp ", im(r"m\leq1-10^{-14}"),
          " (rim only; observation points stay outside the magnet)")],
    [Cell(im(r"\rho=0"), " in the ", im(r"1/\rho"), " factor of ", im(r"H_\rho")),
     Cell("guard ", im(r"\max(\rho,10^{-30})"), "; then force ",
          im(r"H_\rho(\rho<10^{-12}R)=0"))],
    [Cell(im(r"0/0"), " in the Cartesian projection on the axis"),
     Cell("mask ", im(r"\rho>\varepsilon"), ", set ",
          im(r"B_x=B_y=0"), " on the axis")],
    [Cell("Integration singularity at the rim"),
     Cell("observation plane lifted by ", im(r"z_0^{\mathrm{add}}"),
          " (default 1 mm) above the magnet")],
]
tbl = Table(issues, colWidths=[6.0*cm, 9.5*cm])
tbl.setStyle(TableStyle([
    ("BACKGROUND",     (0,0), (-1,0), colors.HexColor("#1a3a6c")),
    ("TEXTCOLOR",      (0,0), (-1,0), colors.white),
    ("VALIGN",         (0,0), (-1,-1), "TOP"),
    ("GRID",           (0,0), (-1,-1), 0.4, colors.HexColor("#cccccc")),
    ("ROWBACKGROUNDS", (0,1), (-1,-1),
                       [colors.white, colors.HexColor("#f4f4f4")]),
    ("LEFTPADDING",    (0,0), (-1,-1), 5),
    ("RIGHTPADDING",   (0,0), (-1,-1), 5),
    ("TOPPADDING",     (0,0), (-1,-1), 4),
    ("BOTTOMPADDING",  (0,0), (-1,-1), 4),
]))
story += [tbl]

# --- 5
story += [
    Paragraph("5. Validation protocol", H1),
    C("Performed each time the analytic core is modified."),
    C("<b>1. On-axis comparison</b> at 200 log-spaced points ",
      im(r"z\in[\,b+10^{-4},\,20L\,]"), ":"),
    C("&bull; numerical (1-D Gauss&ndash;Legendre) vs exact closed form &mdash; "
      "must agree to relative ", im(r"10^{-9}"), " over the full range;"),
    C("&bull; numerical vs dipole far-field &mdash; agree only for ",
      im(r"z\gg L"), "; the log-log error plot must asymptote to a ",
      im(r"1/r^{4}"), "-like decay of the residual."),
    C("<b>2. Divergence check</b> on the 2-D grid: trace of the gradient "
      "tensor"),
    block(r"\nabla\!\cdot\!\mathbf{B}\;=\;"
          r"\partial_x B_x+\partial_y B_y+\partial_z B_z"),
    C("plotted alongside the diagonal components. In source-free conditions "
      "the resulting map should be numerical noise (orders of magnitude "
      "below the component magnitudes)."),
    C("<b>3. Cross-check between programs</b>: at a small set of grid "
      "points, the brute-force script and the analytic script must agree to ",
      im(r"10^{-8}"), " relative tolerance."),
]

# --- 6 : visualization-protocol table
story += [Paragraph("6. Visualization protocol", H1)]
viz = [
    [Cell("<b>Figure</b>"), Cell("<b>Content</b>")],
    [Cell("1"),
     Cell("On-axis ", im(r"B_z"), " and its relative error vs both "
          "references (log&ndash;log, 2 panels)")],
    [Cell("2"),
     Cell(im(r"B_x"), ", ", im(r"B_y"),
          " as <i>pcolor</i> maps, ", im(r"B_z"),
          " as <i>surf</i>, common color limits, axes in cm")],
    [Cell("3"),
     Cell(im(r"B_z"), " color map with quiver of normalized in-plane ",
          im(r"(B_x,B_y)"), " and magnet outline")],
    [Cell("4"),
     Cell("Five gradient components plus ", im(r"\nabla\!\cdot\!\mathbf{B}"),
          " (sanity check) on a 2", im(r"\times"), "3 tiled layout")],
    [Cell("5"),
     Cell("Three components of ", im(r"\nabla|\mathbf{B}|"),
          " obtained from the gradient tensor and ", im(r"\mathbf{B}"),
          " itself")],
]
tbl2 = Table(viz, colWidths=[1.6*cm, 13.9*cm])
tbl2.setStyle(TableStyle([
    ("BACKGROUND",     (0,0), (-1,0), colors.HexColor("#1a3a6c")),
    ("TEXTCOLOR",      (0,0), (-1,0), colors.white),
    ("ALIGN",          (0,0), (0,-1), "CENTER"),
    ("VALIGN",         (0,0), (-1,-1), "TOP"),
    ("GRID",           (0,0), (-1,-1), 0.4, colors.HexColor("#cccccc")),
    ("ROWBACKGROUNDS", (0,1), (-1,-1),
                       [colors.white, colors.HexColor("#f4f4f4")]),
    ("LEFTPADDING",    (0,0), (-1,-1), 5),
    ("RIGHTPADDING",   (0,0), (-1,-1), 5),
    ("TOPPADDING",     (0,0), (-1,-1), 4),
    ("BOTTOMPADDING",  (0,0), (-1,-1), 4),
]))
story += [tbl2]
story += [C("All plots use docked windows, white background, "
            "<i>jet</i> colormap for the field components, default "
            "colormap for the gradients.")]

# --- 7
story += [
    Paragraph("7. Project parameters (current run)", H1),
    Preformatted(
        "R     = 35 mm        cylinder radius\n"
        "L     = 60 mm        cylinder length\n"
        "Br    = 1.0 T        remanent flux density\n"
        "z0add = 1 mm         distance above the top face for the XY map\n"
        "N     = 80           XY grid resolution (analytic version)\n"
        "Nq    = 96           Gauss-Legendre nodes on [0, R]\n"
        "hd    = 1e-6 m       central-difference step for the gradient tensor",
        CODE),
]

# --- 8
story += [
    Paragraph("8. Reproducibility", H1),
    C("Both scripts are self-contained; no toolboxes beyond base MATLAB are "
      "required (<font face='Courier'>ellipke</font> ships with the core "
      "product). Running each script from a clean session reproduces all "
      "figures. The brute-force script is kept as the slow but transparent "
      "reference implementation; the analytic script is the one used for "
      "production runs."),
]

doc.build(story)
print(f"Wrote {OUT_PDF}")
