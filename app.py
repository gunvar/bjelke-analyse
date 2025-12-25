import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# Konstanter for HEB200 og Stål
E = 210e9  # Pa (N/m2)
I_HEB200 = 5696 * 10**-8  # m4 (Moment of inertia)
H_HEB200 = 0.2  # m (Høyde på bjelken)

st.set_page_config(page_title="Bjelkeanalyse v1.0", layout="wide")

st.title("🏗️ Visuell Bjelkeanalyse")
st.markdown("""
Dette verktøyet beregner nedbøyning av en **HEB200** stålbjelke med punktlast i midten. 
*Utviklet som prototype for digital effektivisering.*
""")

# --- SIDEBAR: Parametere ---
st.sidebar.header("Inndata")
L = st.sidebar.slider("Bjelkelengde (m)", 1.0, 15.0, 5.0, 0.5)
P_kn = st.sidebar.number_input("Punktlast i midten (kN)", 0.0, 500.0, 50.0)
opplager = st.sidebar.selectbox("Opplagringstype", ["Fritt opplagret", "Fast innspent"])

P = P_kn * 1000  # Konverter til Newton
x = np.linspace(0, L, 200)

# --- BEREGNINGER ---
def beregn_nedboyning(x, L, P, E, I, type):
    v = np.zeros_like(x)
    # Formler for nedbøyning v(x) for 0 <= x <= L/2 (Symmetrisk)
    mid = L / 2
    for i, xi in enumerate(x):
        temp_x = xi if xi <= mid else L - xi
        if type == "Fritt opplagret":
            # v = (P*x)/(48*E*I) * (3L^2 - 4x^2)
            v[i] = (P * temp_x) / (48 * E * I) * (3 * L**2 - 4 * temp_x**2)
        else: # Fast innspent
            # v = (P*x^2)/(48*E*I) * (3L - 4x)
            v[i] = (P * temp_x**2) / (48 * E * I) * (3 * L - 4 * temp_x)
    return v

v = beregn_nedboyning(x, L, P, E, I_HEB200, opplager)
max_v = np.max(v) * 1000 # mm

# --- VISUALISERING ---
col1, col2 = st.columns([2, 1])

with col1:
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(x, -v * 1000, label="Nedbøyning (mm)", color="#004b87", lw=3) # Norconsult-ish blå
    ax.axhline(0, color='black', lw=1)
    
    # Tegn opplagere
    if opplager == "Fritt opplagret":
        ax.plot([0, L], [0, 0], 'r^', markersize=15, label="Opplager")
    else:
        ax.axvline(0, color='red', lw=5, label="Innspent")
        ax.axvline(L, color='red', lw=5)

    ax.set_xlabel("Lengde på bjelke (m)")
    ax.set_ylabel("Nedbøyning (mm)")
    ax.set_title(f"Profil: HEB200 | Last: {P_kn} kN")
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend()
    st.pyplot(fig)

with col2:
    st.subheader("Resultater")
    st.metric("Maks nedbøyning", f"{max_v:.2f} mm")
    st.info(f"""
    **Tekniske data (HEB200):**
    - E-modul: {E/1e9:.0f} GPa
    - Moment av treghet (I): {I_HEB200*1e8:.0f} cm⁴
    - Egenvekt er ikke medregnet i denne versjonen.
    """)

st.divider()
st.caption("Fremtidige utvidelser: Flere profiler (HEA, Rør), distribuerte laster og utkragere.")
