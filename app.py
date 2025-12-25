import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# --- STYLING OG KONFIGURASJON ---
st.set_page_config(page_title="Bjelkeanalyse Pro v2.0", layout="wide")

# HEX-kode for Streamlit Dark Theme bakgrunn
ST_DARK_BG = "#0e1117"

# Bjelkedata (HEB200)
BEAM_DATA = {
    "HEB200": {"h": 200, "b": 200, "tw": 9, "tf": 15, "A": 78.1, "I": 5696e-8, "E": 210e9}
}

st.title("🏗️ Visuell Bjelkeanalyse Pro")
st.markdown("Dette er verktøyet beregner opplagerkrefter, nedbøyning, moment- og skjærekrefter for en bjelke. Du vil ha anledning til å endre ulike parametere for å se hvordan det gir utslag for bjelken vår!")

# --- SIDEBAR: KONFIGURASJON ---
st.sidebar.header("1. Bjelke & Materiale")
beam_type = st.sidebar.selectbox("Velg bjelketype", list(BEAM_DATA.keys()))
L = st.sidebar.slider("Bjelkelengde (L) [m]", 1.0, 20.0, 8.0, 0.5)

st.sidebar.header("2. Opplagring")
l_supp = st.sidebar.selectbox("Venstre side", ["Fritt opplager", "Fast innspent"])
r_supp = st.sidebar.selectbox("Høyre side", ["Fritt opplager", "Fast innspent"])

st.sidebar.header("3. Belastning")
# Punktlast
use_p = st.sidebar.toggle("Punktlast (P)", value=True)
if use_p:
    P_kn = st.sidebar.number_input("Størrelse (P) [kN]", 0.0, 1000.0, 50.0)
    x_p = st.sidebar.slider("Posisjon [m]", 0.0, L, L/2)
else:
    P_kn, x_p = 0, 0

# Fordelt last
use_q = st.sidebar.toggle("Jevnt fordelt last (q)", value=False)
if use_q:
    q_kn = st.sidebar.number_input("Størrelse (q) [kN/m]", 0.0, 500.0, 10.0)
    q_range = st.sidebar.slider("Utbredelse (fra - til) [m]", 0.0, L, (0.0, L))
else:
    q_kn, q_range = 0, (0, 0)

# --- BEREGNINGER (NUMERISK INTEGRASJON) ---
x = np.linspace(0, L, 500)
E, I = BEAM_DATA[beam_type]["E"], BEAM_DATA[beam_type]["I"]

def get_statics(x_arr, L, P, xp, q, q_r, l_s, r_s):
    # Enkel numerisk løsning ved bruk av superposisjon og diskretisering
    # For å håndtere statisk ubestemte bjelker (fast innspenning),
    # bruker vi "Force Method" forenklet for denne appen.
    
    # Diskretiser lasten
    load = np.zeros_like(x_arr)
    # Legg til punktlast (nærmeste punkt)
    idx_p = (np.abs(x_arr - xp)).argmin()
    load[idx_p] += P / (L/len(x_arr))
    # Legg til fordelt last
    mask_q = (x_arr >= q_r[0]) & (x_arr <= q_r[1])
    load[mask_q] += q
    
    # Skjærkraft V (integralet av last)
    # Her bruker vi en forenklet statisk løsning for reaksjonskrefter
    # (Dette er en forenkling for demo - i en full versjon brukes matrisemetoden)
    # Vi beregner her reaksjonskrefter R1 og R2 basert på momentvekt
    # For punktlast:
    a, b = xp, L - xp
    if l_s == "Fritt opplager" and r_s == "Fritt opplager":
        R1 = P*(b/L) + q*(q_r[1]-q_r[0])*((L-(q_r[0]+q_r[1])/2)/L)
        R2 = (P + q*(q_r[1]-q_r[0])) - R1
    elif l_s == "Fast innspent" and r_s == "Fast innspent":
        R1 = P*b**2*(3*a+b)/L**3 + (q*(q_r[1]-q_r[0]))/2 # Forenklet fordelt
        R2 = (P + q*(q_r[1]-q_r[0])) - R1
    else: # Propped cantilever
        R1 = P*b**2*(3*L-b)/(2*L**3) + (q*(q_r[1]-q_r[0]))*0.6 # Forenklet
        R2 = (P + q*(q_r[1]-q_r[0])) - R1

    # Generer V, M og y diagrammer numerisk
    V = np.zeros_like(x_arr)
    M = np.zeros_like(x_arr)
    y = np.zeros_like(x_arr)
    
    # V(x) = R1 - integral(load)
    V = R1 - np.cumsum(load) * (L/len(x_arr))
    # M(x) = integral(V)
    M = np.cumsum(V) * (L/len(x_arr))
    # y(x) = integral(integral(M)) / EI
    y = np.cumsum(np.cumsum(M)) * (L/len(x_arr))**2 / (E*I)
    y = y - np.linspace(y[0], y[-1], len(y)) # Korreksjon for randbetingelser

    return V, M, y, R1, R2

V, M, y, R1, R2 = get_statics(x, L, P_kn*1000, x_p, q_kn*1000, q_range, l_supp, r_supp)

# --- VISUALISERING ---

def setup_plot():
    fig, ax = plt.subplots(figsize=(6, 1.8))
    fig.patch.set_facecolor(ST_DARK_BG)
    ax.set_facecolor(ST_DARK_BG)
    ax.spines['bottom'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.tick_params(colors='white', labelsize=8)
    ax.grid(True, alpha=0.1, linestyle='--')
    return fig, ax

# 1. Bjelkeskjema
st.subheader("Bjelkeskjema")
fig_s, ax_s = setup_plot()
ax_s.plot([0, L], [0, 0], 'white', lw=3)
# Tegn Opplagere (Ingeniør-stil)
if l_supp == "Fast innspent":
    ax_s.add_patch(patches.Rectangle((-0.2, -0.5), 0.2, 1.0, color='red', alpha=0.7))
    for i in np.linspace(-0.5, 0.4, 5): ax_s.plot([-0.2, 0], [i, i+0.1], 'red', lw=1)
else:
    ax_s.plot([0], [0], 'r^', markersize=12)
if r_supp == "Fast innspent":
    ax_s.add_patch(patches.Rectangle((L, -0.5), 0.2, 1.0, color='red', alpha=0.7))
    for i in np.linspace(-0.5, 0.4, 5): ax_s.plot([L, L+0.2], [i, i+0.1], 'red', lw=1)
else:
    ax_s.plot([L], [0], 'r^', markersize=12)
# Tegn Laster
if use_p: ax_s.annotate('', xy=(x_p, 0), xytext=(x_p, 0.8), arrowprops=dict(facecolor='cyan', shrink=0.05, width=2))
if use_q: ax_s.add_patch(patches.Rectangle((q_range[0], 0), q_range[1]-q_range[0], 0.4, color='orange', alpha=0.3))
# Reaksjonskrefter
ax_s.text(0, -0.8, f"F1: {R1/1000:.1f}kN ↑", color='#4CAF50', fontsize=9, fontweight='bold')
ax_s.text(L, -0.8, f"F2: {R2/1000:.1f}kN ↑", color='#4CAF50', fontsize=9, fontweight='bold', ha='right')
ax_s.set_ylim(-1.2, 1.2)
ax_s.axis('off')
st.pyplot(fig_s)



# 2. Diagrammer (Kompakte)
col1, col2, col3 = st.columns(3)

with col1:
    st.caption("Skjærkraft V [kN]")
    f, a = setup_plot()
    a.fill_between(x, V/1000, color='orange', alpha=0.2)
    a.plot(x, V/1000, color='orange')
    st.pyplot(f)
    st.metric("V_max", f"{np.max(np.abs(V/1000)):.1f} kN")

with col2:
    st.caption("Moment M [kNm]")
    f, a = setup_plot()
    a.fill_between(x, M/1000, color='red', alpha=0.2)
    a.plot(x, M/1000, color='red')
    st.pyplot(f)
    st.metric("M_max", f"{np.max(np.abs(M/1000)):.1f} kNm")

with col3:
    st.caption("Nedbøyning y [mm]")
    f, a = setup_plot()
    a.fill_between(x, -y*1000, color='cyan', alpha=0.2)
    a.plot(x, -y*1000, color='cyan')
    st.pyplot(f)
    st.metric("y_max", f"{np.max(np.abs(y*1000)):.2f} mm")

# --- TEKNISK DATA (BUNN) ---
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
    # Kompakt tverrsnittstegning
    fig_cs, ax_cs = plt.subplots(figsize=(1.5, 1.5))
    fig_cs.patch.set_alpha(0)
    h, b, tw, tf = d['h'], d['b'], d['tw'], d['tf']
    ax_cs.add_patch(patches.Rectangle((-b/2, h/2-tf), b, tf, color='grey'))
    ax_cs.add_patch(patches.Rectangle((-b/2, -h/2), b, tf, color='grey'))
    ax_cs.add_patch(patches.Rectangle((-tw/2, -h/2), tw, h, color='grey'))
    ax_cs.set_xlim(-150, 150); ax_cs.set_ylim(-150, 150); ax_cs.axis('off')
    st.pyplot(fig_cs)
