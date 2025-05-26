# -*- coding: utf-8 -*-
"""
Streamlit app para procesar medidas .s2p con time-gating,
generar 5 figuras iniciales, calcular Axial Ratio y fase,
y exportar resultados.
"""
import streamlit as st
import skrf as rf
import matplotlib.pyplot as plt
import numpy as np
from datetime import date
import io
import tempfile


# Constantes
CLUZ = 299_792_458       # velocidad de la luz (m/s)
R2A  = 180 / np.pi       # radianes a grados

st.title("Procesamiento de Medidas de Red")

# 1) Uploaders para los cuatro archivos .s2p
vv_file  = st.file_uploader("Selecciona Muestra  Antenas Paralelas Rco/Tco/Tyy/VV ", type="s2p")
hv_file  = st.file_uploader("Selecciona Muestra  Antenas Cruzadas  Rcr/Tcr/Tyy/HV ", type="s2p")
vac_file = st.file_uploader("Selecciona Vacío Antenas Paralelas ", type="s2p")
pec_file = st.file_uploader("Selecciona PEC", type="s2p")

# 2) Parámetros de gating y geometría
L_mm      = st.number_input("Espesor de la muestra (mm)", value=7.35)
centerTco = st.number_input("Center Tco (ns)", value=5.2)
spanTco   = st.number_input("Span Tco (ns)",   value=2.0)
centerTcr = st.number_input("Center Tcr (ns)", value=5.25)
spanTcr   = st.number_input("Span Tcr (ns)",   value=2.0)
centerVAC = st.number_input("Center VAC (ns)", value=5.2)
spanVAC   = st.number_input("Span VAC (ns)",   value=2.0)
centerR   = st.number_input("Center R (ns)",   value=6.7)
spanR     = st.number_input("Span R (ns)",     value=2.0)
centerPEC = st.number_input("Center PEC (ns)", value=6.7)
spanPEC   = st.number_input("Span PEC (ns)",   value=2.0)

if st.button("Procesar"):
    # 3) Validación de uploads
    if not (vv_file and hv_file and vac_file and pec_file):
        st.error("Por favor sube los 4 archivos .s2p antes de procesar.")
        st.stop()

    # 4) Escribir cada upload a un fichero temporal y cargar con skrf
    with tempfile.NamedTemporaryFile(suffix=".s2p", delete=False) as tmpvv:
        tmpvv.write(vv_file.getbuffer())
        tmpvv.flush()
        MUESVV = rf.Network(tmpvv.name)
    
    with tempfile.NamedTemporaryFile(suffix=".s2p", delete=False) as tmphv:
        tmphv.write(hv_file.getbuffer())
        tmphv.flush()
        MUESHV = rf.Network(tmphv.name)
    
    with tempfile.NamedTemporaryFile(suffix=".s2p", delete=False) as tmpvac:
        tmpvac.write(vac_file.getbuffer())
        tmpvac.flush()
        vac    = rf.Network(tmpvac.name)
    
    with tempfile.NamedTemporaryFile(suffix=".s2p", delete=False) as tmppec:
        tmppec.write(pec_file.getbuffer())
        tmppec.flush()
        pec    = rf.Network(tmppec.name)
 
    # 5) Time-gating
    MUESTVV_gated = MUESVV.s21.time_gate(
        center=centerTco, span=spanTco, t_unit='ns',
        window=('kaiser',8), method='rfft'
    )
    MUESTHV_gated = MUESHV.s21.time_gate(
        center=centerTcr, span=spanTcr, t_unit='ns',
        window=('kaiser',8), method='rfft'
    )
    MUESRVV_gated = MUESVV.s11.time_gate(
        center=centerR, span=spanR, t_unit='ns',
        window=('kaiser',8), method='rfft'
    )
    vac_gated = vac.s21.time_gate(
        center=centerVAC, span=spanVAC, t_unit='ns',
        window=('kaiser',8), method='rfft'
    )
    pec_gated = pec.s11.time_gate(
        center=centerPEC, span=spanPEC, t_unit='ns',
        window=('kaiser',8), method='rfft'
    )

    # 6) Figuras 1–5

    # Figura 1: MUESVV S21
    fig1, ax1 = plt.subplots(2,1, figsize=(10,6))
    MUESVV.s21.plot_s_time_db(label='raw S21', ax=ax1[0])
    MUESTVV_gated.plot_s_time_db(label='gated S21', ax=ax1[0])
    ax1[0].set_title('MUESVV S21 Tiempo')
    ax1[0].set_xlim(0,10)
    ax1[0].set_ylim(-80,0)
    ax1[0].axvline(centerTco, color='black', linestyle='--')
    ax1[0].text(centerTco+0.2, -80, f'Center={centerTco}', rotation=90)
    ax1[0].legend()
    MUESVV.s21.plot_s_mag(label='raw S21', ax=ax1[1])
    MUESTVV_gated.plot_s_mag(label='gated S21', ax=ax1[1])
    ax1[1].set_title('MUESVV S21 Frecuencia')
    ax1[1].set_ylim(0,1)
    ax1[1].legend()
    plt.tight_layout()
    st.pyplot(fig1)

    # Figura 2: MUESHV S21
    fig2, ax2 = plt.subplots(2,1, figsize=(10,6))
    MUESHV.s21.plot_s_time_db(label='raw S21', ax=ax2[0])
    MUESTHV_gated.plot_s_time_db(label='gated S21', ax=ax2[0])
    ax2[0].set_title('MUESHV S21 Tiempo')
    ax2[0].set_xlim(0,10)
    ax2[0].set_ylim(-80,0)
    ax2[0].axvline(centerTcr, color='black', linestyle='--')
    ax2[0].text(centerTcr+0.2, -80, f'Center={centerTcr}', rotation=90)
    ax2[0].legend()
    MUESHV.s21.plot_s_mag(label='raw S21', ax=ax2[1])
    MUESTHV_gated.plot_s_mag(label='gated S21', ax=ax2[1])
    ax2[1].set_title('MUESHV S21 Frecuencia')
    ax2[1].set_ylim(0,1)
    ax2[1].legend()
    plt.tight_layout()
    st.pyplot(fig2)

    # Figura 3: MUESVV S11
    fig3, ax3 = plt.subplots(2,1, figsize=(10,6))
    MUESVV.s11.plot_s_time_db(label='raw S11', ax=ax3[0])
    MUESRVV_gated.plot_s_time_db(label='gated S11', ax=ax3[0])
    ax3[0].set_title('MUESVV S11 Tiempo')
    ax3[0].set_xlim(0,10)
    ax3[0].set_ylim(-80,0)
    ax3[0].axvline(centerR, color='black', linestyle='--')
    ax3[0].text(centerR+0.2, -80, f'Center={centerR}', rotation=90)
    ax3[0].legend()
    MUESVV.s11.plot_s_mag(label='raw S11', ax=ax3[1])
    MUESRVV_gated.plot_s_mag(label='gated S11', ax=ax3[1])
    ax3[1].set_title('MUESVV S11 Frecuencia')
    ax3[1].set_ylim(0,1)
    ax3[1].legend()
    plt.tight_layout()
    st.pyplot(fig3)

    # Figura 4: VAC S21
    fig4, ax4 = plt.subplots(2,1, figsize=(10,6))
    vac.s21.plot_s_time_db(label='raw S21', ax=ax4[0])
    vac_gated.plot_s_time_db(label='gated S21', ax=ax4[0])
    ax4[0].set_title('VAC S21 Tiempo')
    ax4[0].set_xlim(0,10)
    ax4[0].set_ylim(-80,0)
    ax4[0].axvline(centerVAC, color='black', linestyle='--')
    ax4[0].text(centerVAC+0.2, -80, f'Center={centerVAC}', rotation=90)
    ax4[0].legend()
    vac.s21.plot_s_mag(label='raw S21', ax=ax4[1])
    vac_gated.plot_s_mag(label='gated S21', ax=ax4[1])
    ax4[1].set_title('VAC S21 Frecuencia')
    ax4[1].set_ylim(0,1)
    ax4[1].legend()
    plt.tight_layout()
    st.pyplot(fig4)

    # Figura 5: PEC S11
    fig5, ax5 = plt.subplots(2,1, figsize=(10,6))
    pec.s11.plot_s_time_db(label='raw S11', ax=ax5[0])
    pec_gated.plot_s_time_db(label='gated S11', ax=ax5[0])
    ax5[0].set_title('PEC S11 Tiempo')
    ax5[0].set_xlim(0,10)
    ax5[0].set_ylim(-80,0)
    ax5[0].axvline(centerPEC, color='black', linestyle='--')
    ax5[0].text(centerPEC+0.2, -80, f'Center={centerPEC}', rotation=90)
    ax5[0].legend()
    pec.s11.plot_s_mag(label='raw S11', ax=ax5[1])
    pec_gated.plot_s_mag(label='gated S11', ax=ax5[1])
    ax5[1].set_title('PEC S11 Frecuencia')
    ax5[1].set_ylim(0,1)
    ax5[1].legend()
    plt.tight_layout()
    st.pyplot(fig5)

    # 7) Cálculo de coeficientes en frecuencia
    f       = MUESVV.frequency.f
    fGHz    = f / 1e9
    k0L     = (L_mm/1000) * 2 * np.pi * f / CLUZ

    TVV = MUESTVV_gated.s[:,0,0] / vac_gated.s[:,0,0] * np.exp(-1j * k0L)
    THV = MUESTHV_gated.s[:,0,0] / vac_gated.s[:,0,0] * np.exp(-1j * k0L)
    RVV = MUESRVV_gated.s[:,0,0] / pec_gated.s[:,0,0] * (-1)

    fig_coef, ax_coef = plt.subplots(
        2,1, figsize=(10,6), sharex=True, constrained_layout=True
    )
    ax_coef[0].plot(fGHz, np.abs(TVV), 'b-', fGHz, np.abs(THV), 'r-', fGHz, np.abs(RVV), 'g-')
    ax_coef[0].legend(['TVV','THV','RVV'], loc='best')
    ax_coef[0].set_ylim(0,1)
    ax_coef[0].set_ylabel('Mag')
    ax_coef[0].grid()
    ax_coef[1].plot(fGHz, (np.unwrap(np.angle(TVV)) - np.unwrap(np.angle(THV))) * R2A, 'r-')
    ax_coef[1].legend(['TVV - THV'], loc='best')
    ax_coef[1].set_xlabel('Frequency (GHz)')
    ax_coef[1].set_ylabel('Fase (°)')
    ax_coef[1].grid()
    st.pyplot(fig_coef)

    # 8) Cálculo de Axial Ratio y diferencia de fase
    phase_diff   = np.unwrap(np.angle(TVV)) - np.unwrap(np.angle(THV))
    a            = np.abs(TVV)**4 + np.abs(THV)**4 + 2*np.abs(TVV)**2*np.abs(THV)**2*np.cos(2*phase_diff)
    num          = np.abs(TVV)**2 + np.abs(THV)**2 + np.sqrt(a)
    den          = np.abs(TVV)**2 + np.abs(THV)**2 - np.sqrt(a)
    AR           = np.sqrt(num/den)
    AR_dB        = 20 * np.log10(AR)
    phase_deg    = phase_diff * R2A

    # 9) Mostrar figura final de AR y fase
    fig_ar, ax_ar = plt.subplots(
        2,1, figsize=(10,6), sharex=True, constrained_layout=True
    )
    ax_ar[0].plot(fGHz, AR_dB)
    ax_ar[0].set_title("Axial Ratio vs Frecuencia")
    ax_ar[0].set_ylabel("AR (dB)")
    ax_ar[0].set_ylim(0,80)
    ax_ar[0].grid()
    ax_ar[1].plot(fGHz, AR_dB)
    ax_ar[1].set_xlabel("Frecuencia (GHz)")
    ax_ar[1].set_ylabel("AR (dB)")
    ax_ar[1].set_ylim(0,10)
    ax_ar[1].grid()
    st.pyplot(fig_ar)

    st.success("Procesamiento completado")

    # 10) Exportar resultados a .dat
    header = (
        f"Producido por app_filtradotiempo02.py  fecha:{date.today()}\n"
        "F_GHz\tTVV_mod\tTVV_ang\tTHV_mod\tTHV_ang\tRVV_mod\tRVV_ang\tAR_dB\tPhase_deg\n"
    )
    data = np.column_stack([
        fGHz,
        np.abs(TVV),   np.angle(TVV),
        np.abs(THV),   np.angle(THV),
        np.abs(RVV),   np.angle(RVV),
        AR_dB,         phase_deg
    ])
    buf = io.StringIO()
    np.savetxt(buf, data, fmt="%.4f", delimiter="\t")
    file_content = header + buf.getvalue()

    st.download_button(
        label="Descargar resultados (.dat)",
        data=file_content,
        file_name=f"resultados_{date.today()}.dat",
        mime="text/plain"
    )

