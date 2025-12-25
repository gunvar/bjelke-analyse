import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from fpdf import FPDF
import base64

# --- KONFIGURASJON ---
st.set_page_config(page_title="Bjelkeanalyse Pro v3.0", layout="wide")
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
    if P > 0:
        idx_p = (np.abs(x_arr - xp)).argmin()
        load[idx_p] += P / (L/len(x_arr))
    if q > 0:
        mask_q = (x_arr >= q_r[0]) & (x_arr <= q_r[1])
        load[mask_q] += q
    
    # Forenklet statikk for reaksjonskrefter (Superposisjon)
    a, b = xp, L - xp
    if l_s == "Fritt opplager" and r_s == "Fritt opplager":
        R1 = P*(b/L) + q*(q_r[1]-q_r[0])*((L-(q_r[0]+q_r[1])/2)/L)
    elif l_s == "Fast innspent" and r_s == "Fast innspent":
        R1 = P*b**2*(3*a+b)/L**3 + (q*(q_r[1]-q_r[0]))/2
    else:
        R1 = P*b**2*(3*L-b)/(2*L**3) + (q*(q_r[1]-q_r[0]))*0.6
    
    R2 = (P + q*(q_r[1]-q_r[0])) - R1
    V = R1 - np.cumsum(load) * (L/len(x_arr))
    M = np.cumsum(V) * (L/len(x_arr))
    # Inverter nedbøyning for intuitiv visning (y_vis = -y_faktisk)
    y_raw = np.cumsum(np.cumsum(M)) * (L/len(x_arr))**2 / (E*I)
    y_corr = y_raw - np.linspace(y_raw[0], y_raw[-1], len(y_raw))
    return V, M, y_corr, R1, R2

V, M, y, R1, R2 = get_statics(x, L, P_kn*1000, x_p, q_kn*1000, q_range, l_supp, r_supp)

# --- FORMEL-LOGIKK ---
def get_formula(l_s, r_s, p_on, q_on, q_r, L_val):
    if q_on and (q_r[0] > 0 or q_r[1] < L_val):
        return "For komplisert å vise (delvis last)"
    
    f_M, f_y = "", ""
    if l_s == "Fritt opplager" and r_s == "Fritt opplager":
        if p_on: f_M += "PL/4"; f_y += "PL^3/48EI"
        if q_on: f_M += " + qL^2/8"; f_y += " + 5qL^4/384EI"
    elif l_s == "Fast innspent" and r_s == "Fast innspent":
        if p_on: f_M += "PL/8"; f_y += "PL^3/192EI"
        if q_on: f_M += " + qL^2/12"; f_y += " + qL^4/384EI"
    else:
        f_M, f_y = "Sjekk tabell (Propped)", "Sjekk tabell (Propped)"
    return f_M, f_y

form_M, form_y = get_formula(l_supp, r_supp, use_p, use_q, q_range, L)

# --- PLOTTING ---
def setup_plot(size=(8, 2)):
    fig, ax = plt.subplots(figsize=size)
    fig.patch.set_facecolor(ST_DARK_BG); ax.set_facecolor(ST_DARK_BG)
    ax.spines['bottom'].set_color('white'); ax.spines['left'].set_color('white')
    ax.tick_params(colors='white', labelsize=8); ax.grid(True, alpha=0.1)
    return fig, ax

# 1. Bjelkeskjema (50% mindre enn før)
st.write("### Bjelkeskjema")
fig_s, ax_s = setup_plot(size=(6, 1.2)) # Mindre høyde og bredde
ax_s.plot([0, L], [0, 0], 'white', lw=3)
if l_supp == "Fast innspent":
    ax_s.add_patch(patches.Rectangle((-0.15, -0.4), 0.15, 0.8, color='red', alpha=0.6))
else: ax_s.plot([0], [0], 'r^', markersize=10)
if r_supp == "Fast innspent":
    ax_s.add_patch(patches.Rectangle((L, -0.4), 0.15, 0.8, color='red', alpha=0.6))
else: ax_s.plot([L], [0], 'r^', markersize=10)
if use_p: ax_s.annotate('', xy=(x_p, 0), xytext=(x_p, 0.6), arrowprops=dict(facecolor='cyan', shrink=0.05, width=1.5))
ax_s.text(0, -0.7, f"F1: {R1/1000:.1f}kN", color='#4CAF50', fontsize=8)
ax_s.text(L, -0.7, f"F2: {R2/1000:.1f}kN", color='#4CAF50', fontsize=8, ha='right')
ax_s.axis('off'); st.pyplot(fig_s)

# 2. Vertikale diagrammer
def plot_row(title, data, color, unit, formula, max_val, invert=False):
    col_graph, col_txt = st.columns([4, 1])
    with col_graph:
        f, a = setup_plot()
        display_data = -data if invert else data
        a.fill_between(x, display_data, color=color, alpha=0.2)
        a.plot(x, display_data, color=color)
        st.pyplot(f)
    with col_txt:
        st.write(f"**{title}**")
        st.latex(formula)
        st.metric("Max", f"{max_val:.2f} {unit}")

st.divider()
plot_row("Nedbøyning (y)", y*1000, "cyan", "mm", f"y_{{max}} = {form_y}", np.max(np.abs(y*1000)), invert=True)
plot_row("Moment (M)", M/1000, "red", "kNm", f"M_{{max}} = {form_M}", np.max(np.abs(M/1000)))
plot_row("Skjærkraft (V)", V/1000, "orange", "kN", "V_{max} = R_{max}", np.max(np.abs(V/1000)))

# --- PDF EKSPORT ---
def create_pdf():
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", 'B', 16)
    pdf.cell(200, 10, "Bjelkeanalyse Rapport", ln=True, align='C')
    pdf.set_font("Arial", '', 12)
    pdf.ln(10)
    pdf.cell(200, 10, f"Bjelketype: {beam_type} | Lengde: {L}m", ln=True)
    pdf.cell(200, 10, f"Maks Nedbøyning: {np.max(np.abs(y*1000)):.2f} mm", ln=True)
    pdf.cell(200, 10, f"Maks Moment: {np.max(np.abs(M/1000)):.2f} kNm", ln=True)
    return pdf.output()

st.sidebar.divider()
if st.sidebar.button("📄 Eksporter til PDF"):
    pdf_content = create_pdf()
    st.sidebar.download_button(label="Last ned PDF", data=pdf_content, file_name="Bjelkeanalyse.pdf", mime="application/pdf")
