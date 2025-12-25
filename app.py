import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# --- KONFIGURASJON OG DATA ---
st.set_page_config(page_title="Bjelkeanalyse Pro", layout="wide")

# Bjelkedata (Kan utvides senere)
BEAM_DATA = {
    "HEB200": {
        "h": 200, "b": 200, "tw": 9, "tf": 15, 
        "A": 78.1, "I": 5696 * 10**-8, "E": 210e9,
        "weight": 61.3
    }
}

# --- INTRODUKSJON ---
st.title("🏗️ Visuell Bjelkeanalyse Pro")
st.markdown("Dette er verktøyet beregner opplagerkrefter, nedbøyning, moment- og skjærekrefter for en bjelke. Du vil ha anledning til å endre ulike parametere for å se hvordan det gir utslag for bjelken vår!")

# --- SIDEBAR (INPUT) ---
st.sidebar.header("Konfigurasjon")
beam_type = st.sidebar.selectbox("Velg bjelketype", list(BEAM_DATA.keys()))
L = st.sidebar.slider("Bjelkelengde (L) [m]", 1.0, 15.0, 6.0, 0.5)
P_kn = st.sidebar.number_input("Punktlast i midten (P) [kN]", 0.0, 500.0, 100.0)

st.sidebar.subheader("Opplagring")
left_supp = st.sidebar.selectbox("Venstre side", ["Fritt opplager", "Fast innspent"])
right_supp = st.sidebar.selectbox("Høyre side", ["Fritt opplager", "Fast innspent"])

# Hent tekniske data
data = BEAM_DATA[beam_type]
E, I, P = data["E"], data["I"], P_kn * 1000
x = np.linspace(0, L, 500)

# --- STATIKK-MOTOR ---
# Beregner basert på kombinasjon (Punktlast i midten a = L/2)
def calculate_statics(x, L, P, E, I, l_type, r_type):
    V = np.zeros_like(x)
    M = np.zeros_like(x)
    y = np.zeros_like(x)
    mid = L/2
    
    if l_type == "Fritt opplager" and r_type == "Fritt opplager":
        R1 = R2 = P/2
        M_max_val = (P*L)/4
        y_max_val = (P*L**3)/(48*E*I)
        f_V, f_M, f_y = "P/2", "PL/4", "PL^3 / 48EI"
        for i, xi in enumerate(x):
            V[i] = R1 if xi < mid else -R2
            M[i] = (P*xi)/2 if xi <= mid else (P*(L-xi))/2
            y[i] = (P*xi*(3*L**2 - 4*xi**2))/(48*E*I) if xi <= mid else (P*(L-xi)*(3*L**2 - 4*(L-xi)**2))/(48*E*I)
            
    elif l_type == "Fast innspent" and r_type == "Fast innspent":
        R1 = R2 = P/2
        M_max_val = (P*L)/8
        y_max_val = (P*L**3)/(192*E*I)
        f_V, f_M, f_y = "P/2", "PL/8", "PL^3 / 192EI"
        for i, xi in enumerate(x):
            V[i] = R1 if xi < mid else -R2
            M[i] = (P/8)*(4*xi - L) if xi <= mid else (P/8)*(3*L - 4*xi)
            y[i] = (P*xi**2 * (3*L - 4*xi))/(48*E*I) if xi <= mid else (P*(L-xi)**2 * (3*L - 4*(L-xi)))/(48*E*I)

    else: # En side fast, en side fritt (Propped cantilever)
        # Antar venstre er fast hvis ulikt, ellers speiles det logisk
        is_left_fixed = (l_type == "Fast innspent")
        R_fixed = (11*P)/16
        R_pinned = (5*P)/16
        R1, R2 = (R_fixed, R_pinned) if is_left_fixed else (R_pinned, R_fixed)
        M_max_val = (3*P*L)/16
        y_max_val = (7*P*L**3)/(768*E*I)
        f_V, f_M, f_y = "11P/16", "3PL/16", "7PL^3 / 768EI"
        for i, xi in enumerate(x):
            pos = xi if is_left_fixed else L-xi
            V[i] = R1 if xi < mid else -R2
            if pos <= mid:
                M_val = (P/16)*(11*pos - 3*L)
                y_val = (P*pos**2 * (9*L - 11*pos))/(96*E*I)
            else:
                M_val = (5*P/16)*(L - pos)
                y_val = (P*(L-pos)*(3*L**2 - 5*(L-pos)**2))/(96*E*I)
            M[i] = M_val
            y[i] = y_val

    return x, V, M, y, R1, R2, f_V, f_M, f_y, M_max_val, y_max_val

x, V, M, y, R1, R2, f_V, f_M, f_y, M_max_val, y_max_val = calculate_statics(x, L, P, E, I, left_supp, right_supp)

# --- 1. SKJEMA-FIGUR (TOPP) ---
fig0, ax0 = plt.subplots(figsize=(10, 2))
ax0.plot([0, L], [0, 0], 'black', lw=4) # Bjelke
# Tegn opplagere
if left_supp == "Fast innspent": ax0.axvline(0, color='red', lw=10)
else: ax0.plot(0, 0, 'r^', markersize=15)
if right_supp == "Fast innspent": ax0.axvline(L, color='red', lw=10)
else: ax0.plot(L, 0, 'r^', markersize=15)
# Last-pil
ax0.annotate('', xy=(L/2, 0), xytext=(L/2, 1), arrowprops=dict(facecolor='blue', shrink=0.05))
ax0.text(L/2, 1.1, f'P = {P_kn}kN', ha='center')
# Opplagerkrefter
ax0.text(0, -0.8, f'F1={R1/1000:.1f}kN ↑', color='green', fontweight='bold')
ax0.text(L, -0.8, f'F2={R2/1000:.1f}kN ↑', color='green', fontweight='bold', ha='right')
ax0.set_ylim(-1.5, 2)
ax0.axis('off')
st.pyplot(fig0)

# --- 2. DIAGRAMMER ---
def plot_diagram(x, values, title, unit, formula, max_val, color):
    col_a, col_b = st.columns([3, 1])
    with col_a:
        fig, ax = plt.subplots(figsize=(8, 2.5))
        ax.plot(x, values, color=color, lw=2)
        ax.fill_between(x, values, color=color, alpha=0.1)
        ax.axhline(0, color='black', lw=1)
        ax.set_ylabel(f"[{unit}]")
        ax.grid(True, alpha=0.3)
        st.pyplot(fig)
    with col_b:
        st.write(f"**{title}**")
        st.latex(f"{formula}")
        st.metric("Maks verdi", f"{abs(max_val):.2f} {unit}")

st.divider()
plot_diagram(x, V/1000, "Skjærkraft (V)", "kN", f"V_{{max}} = {f_V}", max(abs(V/1000)), "orange")
plot_diagram(x, M/1000, "Bøyemoment (M)", "kNm", f"M_{{max}} = {f_M}", M_max_val/1000, "red")
plot_diagram(x, -y*1000, "Nedbøyning (y)", "mm", f"y_{{max}} = {f_y}", y_max_val*1000, "blue")

# --- 3. TEKNISK DATA (BUNN) ---
st.divider()
with st.container():
    st.subheader(f"Tekniske data: {beam_type}")
    c1, c2, c3 = st.columns([1, 1, 2])
    with c1:
        st.write(f"**Materiale:** Stål")
        st.write(f"**E-modul:** {data['E']/1e9} GPa")
        st.write(f"**Areal (A):** {data['A']} cm²")
        st.write(f"**Treghetsmoment (I):** {data['I']*1e8:.0f} cm⁴")
    with c2:
        st.write(f"**Høyde (h):** {data['h']} mm")
        st.write(f"**Bredde (b):** {data['b']} mm")
        st.write(f"**Steg (tw):** {data['tw']} mm")
        st.write(f"**Flens (tf):** {data['tf']} mm")
    with c3:
        # Enkel tegning av HEB-tverrsnitt
        fig_cs, ax_cs = plt.subplots(figsize=(2, 2))
        h, b, tw, tf = data['h'], data['b'], data['tw'], data['tf']
        # Flenser og steg
        ax_cs.add_patch(plt.Rectangle((-b/2, h/2-tf), b, tf, color='grey')) # Topp
        ax_cs.add_patch(plt.Rectangle((-b/2, -h/2), b, tf, color='grey'))   # Bunn
        ax_cs.add_patch(plt.Rectangle((-tw/2, -h/2), tw, h, color='grey'))  # Steg
        ax_cs.set_xlim(-150, 150); ax_cs.set_ylim(-150, 150)
        ax_cs.axis('off')
        st.pyplot(fig_cs)
