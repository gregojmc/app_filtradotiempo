# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 12:41:49 2025
@author: gregomc

Para lanzar:
    streamlit run app_filtradotiempo02.py
"""
import streamlit as st
import skrf as rf
import matplotlib.pyplot as plt
import os
import numpy as np
from datetime import date
import io

ThisCode = 'app_filtradotiempo02.py'
cluz = 299792458  # velocidad de la luz (m/s)
r2a = 180/np.pi   # radianes a grados

st.title("Procesamiento de Medidas de Red")

# Path default o seleccionado
default_path = ""
path = st.text_input("Path (directorio completo)", value=default_path)
L = st.number_input("Espesor de la muestra, (mm)", value=7.35)

# Nombre de la muestra
muestra = st.text_input("Nombre de la muestra")

# Centros y span
col1, col2 = st.columns(2)
with col1:
    centerTco = st.number_input("Center Tco", value=5.2)
    centerTcr = st.number_input("Center Tcr", value=5.25)
    centerVAC = st.number_input("Center VAC", value=5.2)
    centerR   = st.number_input("Center R",   value=6.7)
    centerPEC = st.number_input("Center PEC", value=6.7)
with col2:
    spanTco = st.number_input("Span Tco", value=2.0)
    spanTcr = st.number_input("Span Tcr", value=2.0)
    spanVAC = st.number_input("Span VAC", value=2.0)
    spanR   = st.number_input("Span R",   value=2.0)
    spanPEC = st.number_input("Span PEC", value=2.0)

if st.button('Procesar'):
    if not muestra:
        st.error("Por favor, ingrese el nombre de la muestra")
    else:
        # 1) Montar rutas de los ficheros
        file_paths = {
            'MUESVV': os.path.join(path, f"{muestra}_P1TxV_P2RxV.s2p"),
            'MUESHV': os.path.join(path, f"{muestra}_P1TxV_P2RxH.s2p"),
            'vac':    os.path.join(path,             "VAC_P1TxV_P2RxV.s2p"),
            'pec':    os.path.join(path,             "PEC_P1TxV_P2RxV.s2p")
        }
        # 2) Comprobar existencia
        missing_paths = [fp for fp in file_paths.values() if not os.path.exists(fp)]
        if missing_paths:
            # Mostrar error con la lista de rutas faltantes
            st.error(f"Archivos no encontrados en: {missing_paths}")
            # Mostrar además qué variable corresponde a cada ruta
            for name, fp in file_paths.items():
                if fp in missing_paths:
                    st.write(f"{name}: {fp}")
            st.stop()
        
        # 3) Mostrar rutas (opcional)
        for name, fp in file_paths.items():
            st.write(f"Ruta {name}: {fp}")

        try:
            with st.spinner('Procesando datos...'):
                # Carga de redes
                MUESVV = rf.Network(file_paths['MUESVV'])
                MUESHV = rf.Network(file_paths['MUESHV'])
                vac    = rf.Network(file_paths['vac'])
                pec    = rf.Network(file_paths['pec'])

                # Time gating
                MUESTVV_gated = MUESVV.s21.time_gate(
                    center=centerTco, span=spanTco, t_unit='ns', window=('kaiser',8), method='rfft'
                )
                MUESTHV_gated = MUESHV.s21.time_gate(
                    center=centerTcr, span=spanTcr, t_unit='ns', window=('kaiser',8), method='rfft'
                )
                MUESRVV_gated = MUESVV.s11.time_gate(
                    center=centerR, span=spanR, t_unit='ns', window=('kaiser',8), method='rfft'
                )
                vac_gated = vac.s21.time_gate(
                    center=centerVAC, span=spanVAC, t_unit='ns', window=('kaiser',8), method='rfft'
                )
                pec_gated = pec.s11.time_gate(
                    center=centerPEC, span=spanPEC, t_unit='ns', window=('kaiser',8), method='rfft'
                )

                # Figuras 1–5
                # Figura 1
                fig1, ax1 = plt.subplots(2,1, figsize=(10,6))
                MUESVV.s21.plot_s_time_db(label='raw S21', ax=ax1[0])
                MUESTVV_gated.plot_s_time_db(label='gated S21', ax=ax1[0])
                ax1[0].set_title('MUESVV S21 Tiempo')
                ax1[0].set_xlim(0,10); ax1[0].set_ylim(-80,0)
                ax1[0].axvline(centerTco, color='black', ls='--')
                ax1[0].text(centerTco+0.2, -80, f'Center={centerTco}', rotation=90)
                ax1[0].legend()
                MUESVV.s21.plot_s_mag(label='raw S21', ax=ax1[1])
                MUESTVV_gated.plot_s_mag(label='gated S21', ax=ax1[1])
                ax1[1].set_title(f'MUESVV S21 Frecuencia - {muestra}')
                ax1[1].set_ylim(0,1); ax1[1].legend()
                plt.tight_layout(); st.pyplot(fig1)

                # Figura 2
                fig2, ax2 = plt.subplots(2,1, figsize=(10,6))
                MUESHV.s21.plot_s_time_db(label='raw S21', ax=ax2[0])
                MUESTHV_gated.plot_s_time_db(label='gated S21', ax=ax2[0])
                ax2[0].set_title('MUESHV S21 Tiempo')
                ax2[0].set_xlim(0,10); ax2[0].set_ylim(-80,0)
                ax2[0].axvline(centerTcr, color='black', ls='--')
                ax2[0].text(centerTcr+0.2, -80, f'Center={centerTcr}', rotation=90)
                ax2[0].legend()
                MUESHV.s21.plot_s_mag(label='raw S21', ax=ax2[1])
                MUESTHV_gated.plot_s_mag(label='gated S21', ax=ax2[1])
                ax2[1].set_title('MUESHV S21 Frecuencia'); ax2[1].set_ylim(0,1); ax2[1].legend()
                plt.tight_layout(); st.pyplot(fig2)

                # Figura 3
                fig3, ax3 = plt.subplots(2,1, figsize=(10,6))
                MUESVV.s11.plot_s_time_db(label='raw S11', ax=ax3[0])
                MUESRVV_gated.plot_s_time_db(label='gated S11', ax=ax3[0])
                ax3[0].set_title('MUESVV S11 Tiempo')
                ax3[0].set_xlim(0,10); ax3[0].set_ylim(-80,0)
                ax3[0].axvline(centerR, color='black', ls='--')
                ax3[0].text(centerR+0.2, -80, f'Center={centerR}', rotation=90)
                ax3[0].legend()
                MUESVV.s11.plot_s_mag(label='raw S11', ax=ax3[1])
                MUESRVV_gated.plot_s_mag(label='gated S11', ax=ax3[1])
                ax3[1].set_title('MUESVV S11 Frecuencia'); ax3[1].set_ylim(0,1); ax3[1].legend()
                plt.tight_layout(); st.pyplot(fig3)

                # Figura 4
                fig4, ax4 = plt.subplots(2,1, figsize=(10,6))
                vac.s21.plot_s_time_db(label='raw S21', ax=ax4[0])
                vac_gated.plot_s_time_db(label='gated S21', ax=ax4[0])
                ax4[0].set_title('VAC S21 Tiempo')
                ax4[0].set_xlim(0,10); ax4[0].set_ylim(-80,0)
                ax4[0].axvline(centerVAC, color='black', ls='--')
                ax4[0].text(centerVAC+0.2, -80, f'Center={centerVAC}', rotation=90)
                ax4[0].legend()
                vac.s21.plot_s_mag(label='raw S21', ax=ax4[1])
                vac_gated.plot_s_mag(label='gated S21', ax=ax4[1])
                ax4[1].set_title('VAC S21 Frecuencia'); ax4[1].set_ylim(0,1); ax4[1].legend()
                plt.tight_layout(); st.pyplot(fig4)

                # Figura 5
                fig5, ax5 = plt.subplots(2,1, figsize=(10,6))
                pec.s11.plot_s_time_db(label='raw S11', ax=ax5[0])
                pec_gated.plot_s_time_db(label='gated S11', ax=ax5[0])
                ax5[0].set_title('PEC S11 Tiempo')
                ax5[0].set_xlim(0,10); ax5[0].set_ylim(-80,0)
                ax5[0].axvline(centerPEC, color='black', ls='--')
                ax5[0].text(centerPEC+0.2, -80, f'Center={centerPEC}', rotation=90)
                ax5[0].legend()
                pec.s11.plot_s_mag(label='raw S11', ax=ax5[1])
                pec_gated.plot_s_mag(label='gated S11', ax=ax5[1])
                ax5[1].set_title('PEC S11 Frecuencia'); ax5[1].set_ylim(0,1); ax5[1].legend()
                plt.tight_layout(); st.pyplot(fig5)

                # Cálculo de R y T
                f    = MUESVV.frequency.f
                fGHz = f/1e9
                k0L  = L * 2*np.pi*f/cluz

                TVV = MUESTVV_gated.s[:,0,0]/vac_gated.s[:,0,0]*np.exp(-1j*k0L)
                THV = MUESTHV_gated.s[:,0,0]/vac_gated.s[:,0,0]*np.exp(-1j*k0L)
                RVV = MUESRVV_gated.s[:,0,0]/pec_gated.s[:,0,0]*(-1)

                fig_coef, ax_coef = plt.subplots(2,1, figsize=(10,6), sharex=True, constrained_layout=True)
                ax_coef[0].plot(fGHz, np.abs(TVV), 'b-', fGHz, np.abs(THV), 'r-', fGHz, np.abs(RVV), 'g-')
                ax_coef[0].legend(['TVV','THV','RVV'], loc='best')
                ax_coef[0].set_ylim(0,1); ax_coef[0].set_ylabel('R & T [mag]'); ax_coef[0].grid()
                ax_coef[1].plot(fGHz, (np.unwrap(np.angle(TVV))-np.unwrap(np.angle(THV)))*r2a, 'r-')
                ax_coef[1].legend(['TVV - THV'], loc='best')
                ax_coef[1].set_xlabel('Frequency (GHz)'); ax_coef[1].set_ylabel('Angle [deg]'); ax_coef[1].grid()
                st.pyplot(fig_coef)

                # --- Cálculo Axial Ratio y Diferencia de Fase ---
                phase_diff = np.unwrap(np.angle(TVV)) - np.unwrap(np.angle(THV))
                a = np.abs(TVV)**4 + np.abs(THV)**4 + 2*np.abs(TVV)**2*np.abs(THV)**2*np.cos(2*phase_diff)
                num = np.abs(TVV)**2 + np.abs(THV)**2 + np.sqrt(a)
                den = np.abs(TVV)**2 + np.abs(THV)**2 - np.sqrt(a)
                AR = np.sqrt(num/den)
                AR_dB = 20 * np.log10(AR)
                phase_diff_deg = phase_diff * r2a

                # Figura final
                fig_ar, ax_ar = plt.subplots(2,1, figsize=(10,6), sharex=True, constrained_layout=True)
                ax_ar[0].plot(fGHz, AR_dB)
                ax_ar[0].set_title("Axial Ratio vs Frecuencia")
                ax_ar[0].set_ylabel("AR (dB)");  ax_ar[0].grid()
                ax_ar[1].plot(fGHz, AR_dB)
                ax_ar[1].set_ylabel("AR (dB)")
                ax_ar[1].grid(True); 
                ax_ar[1].set_ylim(0,10)
                st.pyplot(fig_ar)

                # Éxito
                st.success("Procesamiento completado")

                # Preparar salida
                cabecero = (
                    f"Producido por: {ThisCode} fecha:{date.today()}\n"
                    f"Medidas: {path} {muestra}\n"
                    f"Center (Tco,Tcr,R,VAC,PEC): "
                    f"{centerTco:.3f}, {centerTcr:.3f}, {centerR:.3f}, "
                    f"{centerVAC:.3f}, {centerPEC:.3f}\n"
                    f"Span (Tco,Tcr,R,VAC,PEC): "
                    f"{spanTco:.3f}, {spanTcr:.3f}, {spanR:.3f}, "
                    f"{spanVAC:.3f}, {spanPEC:.3f}\n"
                    "F_GHz\tTVV_mod\tTVV_ang\tTHV_mod\tTHV_ang\tRVV_mod\tRVV_ang\tAR_dB\tPhase_diff_deg"
                )

                datos_salida = np.column_stack((
                    fGHz,
                    np.abs(TVV),        np.angle(TVV),
                    np.abs(THV),        np.angle(THV),
                    np.abs(RVV),        np.angle(RVV),
                    AR_dB,              phase_diff_deg
                ))

                buf = io.StringIO()
                np.savetxt(buf, datos_salida, fmt=' %.4f')
                contenido = buf.getvalue()
                archivo_completo = cabecero + "\n" + contenido

                output_filename = st.text_input("Nombre del fichero de salida", value=f"{muestra}_resultados.dat")
                st.download_button(
                    label="Descargar resultados",
                    data=archivo_completo,
                    file_name=output_filename,
                    mime="text/plain"
                )

                # Guardar en disco
                fileout = os.path.join(path, output_filename)
                try:
                    with open(fileout, 'w') as f:
                        f.write(archivo_completo)
                    st.success(f"Archivo también guardado en: {fileout}")
                except Exception as e:
                    st.warning(f"No se pudo guardar en disco: {e}")

        except Exception as e:
            st.error(f"Error en el procesamiento: {e}")
            import traceback
            st.code(traceback.format_exc())
