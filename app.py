import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from fpdf import FPDF

# --- KONFIGURASJON ---
st.set_page_config(page_title="Bjelkeanalyse Pro v3.1", layout="wide")
ST_DARK_BG = "#0e1117"

BEAM_DATA = {
    "HEB200": {"h": 200, "b": 200, "tw": 9, "tf": 15, "A": 78.1, "I": 5696e-8, "E": 210e9}
}

st.title("🏗️ Visuell Bjelkeanalyse Pro")
st.markdown("Dette er verktøyet beregner opplagerkrefter, nedbøyning, moment- og skjærekrefter for en bjelke.")

# --- SIDEBAR ---
st.sidebar.header("Parametere")
beam_type = st.sidebar.selectbox("Velg bjelketype", list(BEAM_DATA.keys()))
L = st.sidebar.slider("Bjelkelengde (L) [m]", 1.0, 20.0, 8.0, 0.5)

st.sidebar.subheader("Opplagring")
l_supp = st.sidebar.selectbox("Venstre side", ["Fritt opplager", "Fast innspent"])
r_supp = st.sidebar.selectbox("Høyre side", ["Fritt opplager", "Fast innspent"])

st.sidebar.subheader("Belastning")
use_p = st.sidebar.toggle("Punktlast (P)", value=True)
P_kn = st.sidebar.number_input("Størrelse (P) [kN]", 0.0, 1000.0, 50.0) if use_p else 0
x_p = st.sidebar.slider("Posisjon [m]", 0.0, L, L/2) if use_p else 0

use_q = st.sidebar.toggle("Jevnt fordelt last (q)", value=False)
q_kn = st.sidebar.number_input("Størrelse (q) [kN/m]", 0.0, 500.0, 10.0) if use_q else 0
q_range = st.sidebar.slider("Utbredelse [m]", 0.0, L, (0.0, L)) if use_q else (0, 0)

# --- BEREGNINGER ---
x = np.linspace(0, L, 500)
E, I = BEAM_DATA[beam_type]["E"], BEAM_DATA[beam_type]["I"]

def get_statics(x_arr, L, P, xp, q, q_r, l_s, r_s):
    load = np.zeros_like(x_arr)
    dx = L/len(x_arr)
    if P > 0:
        idx_p = (np.abs(x_arr - xp)).argmin()
        load[idx_p] += P / dx
    if q > 0:
        mask_q = (x_arr >= q_r[0]) & (x_arr <= q_r[1])
        load[mask_q] += q
    
    a, b = xp, L - xp
    # Beregning av reaksjonskrefter basert på kombinasjon
    if l_s == "Fritt opplager" and r_s == "Fritt opplager":
        R1 = P*(b/L) + (q*(q_r[1]-q_r[0])*(L-(q_r[0]+q_r[1])/2))/L
    elif l_s == "Fast innspent" and r_s == "Fast innspent":
        R1 = P*b**2*(3*a+b)/L**3 + (q*(q_r[1]-q_r[0]))/2 # Forenklet for q
    elif l_s == "Fast innspent" and r_s == "Fritt opplager":
        R1 = (P*b/(2*L**3))*(3*L**2-b**2) + (q*(q_r[1]-q_r[0]))*0.625
    else: # Fritt venstre, fast høyre
        R2 = (P*a/(2*L**3))*(3*L**2-a**2) + (q*(q_r[1]-q_r[0]))*0.625
        R1 = (P + q*(q_r[1]-q_r[0])) - R2
    
    R2 = (P + q*(q_r[1]-q_r[0])) - R1
    V = R1 - np.cumsum(load) * dx
    M = np.cumsum(V) * dx
    # Nedbøyning (E*I*y'' = M). Vi integrerer to ganger og korrigerer for opplager
    y_raw = np.cumsum(np.cumsum(M)) * dx**2 / (E*I)
    y_corr = y_raw - np.linspace(y_raw[0], y_raw[-1], len(y_raw))
    return V, M, y_corr, R1, R2

V, M, y, R1, R2 = get_statics(x, L, P_kn*1000, x_p, q_kn*1000, q_range, l_supp, r_supp)

# --- FORMEL-LOGIKK ---
def get_formula(l_s, r_s, p_on, q_on, q_r, L_v):
    if q_on and (q_r[0] > 0 or q_r[1] < L_v): return "Delvis last: Se diagram", "Delvis last: Se diagram"
    
    # Punktlast i midten (L/2) formler
    if l_s == "Fritt opplager" and r_s == "Fritt opplager":
        fM = "PL/4" if p_on else ""; fy = "PL^3/48EI" if p_on else ""
        if q_on: fM += " + qL^2/8"; fy += " + 5qL^4/384EI"
    elif l_s == "Fast innspent" and r_s == "Fast innspent":
        fM = "PL/8" if p_on else ""; fy = "PL^3/192EI" if p_on else ""
        if q_on: fM += " + qL^2/12"; fy += " + qL^4/384EI"
    else: # Propped (Fast-Fritt)
        fM = "3PL/16" if p_on else ""; fy = "7PL^3/768EI" if p_on else ""
        if q_on: fM += " + qL^2/8"; fy += " + qL^4/185EI"
    return fM, fy

form_M, form_y = get_formula(l_supp, r_supp, use_p, use_q, q_range, L)

# --- PLOTTING ---
def setup_plot(size=(8, 2)):
    fig, ax = plt.subplots(figsize=size)
    fig.patch.set_facecolor(ST_DARK_BG); ax.set_facecolor(ST_DARK_BG)
    ax.spines['bottom'].set_color('white'); ax.spines['left'].set_color('white')
    ax.tick_params(colors='white', labelsize=8); ax.grid(True, alpha=0.1)
    return fig, ax

st.write("### Bjelkeskjema")
fig_s, ax_s = setup_plot(size=(4, 0.8)) # Enda mindre figur
ax_s.plot([0, L], [0, 0], 'white', lw=3)
if l_supp == "Fast innspent": ax_s.add_patch(patches.Rectangle((-0.1, -0.3), 0.1, 0.6, color='red', alpha=0.6))
else: ax_s.plot([0], [0], 'r^', markersize=8)
if r_supp == "Fast innspent": ax_s.add_patch(patches.Rectangle((L, -0.3), 0.1, 0.6, color='red', alpha=0.6))
else: ax_s.plot([L], [0], 'r^', markersize=8)
if use_p: ax_s.annotate('', xy=(x_p, 0), xytext=(x_p, 0.4), arrowprops=dict(facecolor='cyan', shrink=0.05, width=1))
ax_s.text(0, -0.6, f"F1: {R1/1000:.1f}kN", color='#4CAF50', fontsize=7)
ax_s.text(L, -0.6, f"F2: {R2/1000:.1f}kN", color='#4CAF50', fontsize=7, ha='right')
ax_s.axis('off'); st.pyplot(fig_s)

def plot_row(title, data, color, unit, formula, max_val, invert=False):
    col_graph, col_txt = st.columns([4, 1])
    with col_graph:
        f, a = setup_plot(size=(8, 1.8))
        plot_data = -data if invert else data
        a.fill_between(x, plot_data, color=color, alpha=0.2); a.plot(x, plot_data, color=color)
        if invert: a.set_ylim(bottom=min(plot_data)*1.1, top=0) # Tvinger grafen under aksen
        st.pyplot(f)
    with col_txt:
        st.write(f"**{title}**"); st.latex(formula)
        st.metric("Maks", f"{max_val:.2f} {unit}")

st.divider()
plot_row("Nedbøyning (y)", y*1000, "cyan", "mm", f"y_{{max}} \approx {form_y}", np.max(np.abs(y*1000)), invert=True)
plot_row("Bøyemoment (M)", M/1000, "red", "kNm", f"M_{{max}} = {form_M}", np.max(np.abs(M/1000)))
plot_row("Skjærkraft (V)", V/1000, "orange", "kN", "V_{max} = R_{max}", np.max(np.abs(V/1000)))

# --- TEKNISK DATA ---
st.divider()
st.subheader(f"Tekniske data: {beam_type}")
d = BEAM_DATA[beam_type]
c1, c2, c3 = st.columns([1, 1, 1.5])
with c1:
    st.write(f"**E-modul:** {d['E']/1e9} GPa")
    st.write(f"**Treghetsmoment (I):** {d['I']*1e8:.0f} cm⁴")
with c2:
    st.write(f"**Høyde (h):** {d['h']} mm")
    st.write(f"**Bredde (b):** {d['b']} mm")
with c3:
    fig_cs, ax_cs = plt.subplots(figsize=(1.2, 1.2))
    fig_cs.patch.set_alpha(0)
    ax_cs.add_patch(patches.Rectangle((-d['b']/2, d['h']/2-d['tf']), d['b'], d['tf'], color='grey'))
    ax_cs.add_patch(patches.Rectangle((-d['b']/2, -d['h']/2), d['b'], d['tf'], color='grey'))
    ax_cs.add_patch(patches.Rectangle((-d['tw']/2, -d['h']/2), d['tw'], d['h'], color='grey'))
    ax_cs.set_xlim(-150, 150); ax_cs.set_ylim(-150, 150); ax_cs.axis('off'); st.pyplot(fig_cs)

# --- PDF EKSPORT ---
def create_pdf_bytes():
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Helvetica", 'B', 16); pdf.cell(200, 10, "Bjelkeanalyse Rapport", ln=True, align='C')
    pdf.set_font("Helvetica", '', 12); pdf.ln(10)
    pdf.cell(200, 10, f"Bjelketype: {beam_type} | Lengde: {L}m", ln=True)
    pdf.cell(200, 10, f"Maks Nedbøyning: {np.max(np.abs(y*1000)):.2f} mm", ln=True)
    pdf.cell(200, 10, f"Maks Moment: {np.max(np.abs(M/1000)):.2f} kNm", ln=True)
    return pdf.output()

st.sidebar.divider()
if st.sidebar.button("📄 Generer PDF"):
    pdf_data = create_pdf_bytes()
    st.sidebar.download_button(label="Last ned PDF", data=pdf_data, file_name="Bjelkeanalyse.pdf", mime="application/pdf")
