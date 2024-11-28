# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 21:26:16 2023

Das Programm ist eine einfache Monte-Carlo-Simulation für Photonen, die einen
Absorber durchdringen und auf ein ideal absorbierenden Detektor 
(100% Ansprechvermögen und direkte Ionisation und perfekte Energieauflösung).

Der betrachtete Raum ist absolut evakuiert. Die Quelle emittiert die Photonen
isotrop in alle Raumrichtungen. Der Absorber ist unendlich in y- und z-Richtung
ausgedehnt.

Zur Vereinfachung werden nur die inkohärente Streuung und der Photoelektrische
Effekt berücksichtigt. Zusätzlich wird die Berechnung bei Photonen 
unter einer Energie von 0.01 MeV abgebrochen und als absorbiert betrachtet.

@author: NR
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys


# =============================================================================
# Methoden
# =============================================================================

def KN(mu,E):
    """
    Berechnet den differentiellen Wirkungsquerschnitt für einen Streuwinkel bei
    einer Energie E
    """
    gamma = E / 0.511
    return re**2/2*1/(1+gamma*(1-mu))**2 * (1 + mu**2 + 
            gamma**2 * (1 - mu)**2/ (1+gamma*(1-mu)))

def Verwerfungsmethode(E):
    """
    die Verwerfungsmethode dient dazu eine Normalverteilte Zufallsvariable (mu)
    an die Wahrscheinlichkeitsdichtefunktion des differentiellen Wirkungs-
    querschnitts (Formel von Klein-Nishina) anzupassen.    
    """
    while True:
        mu = np.random.uniform(-1, 1) # Zufallszahl für mu (cos(Streuwinkel))
        y = np.random.uniform(0, 1e-29) #obere Grenze = Integral KN für 10keV
        kn_mu = KN(mu,E)
        if y >= kn_mu:          #Filtert die Werte für mu, wenn diese    
            return mu           #unter differentiellen KN liegen

def clear_line():
    sys.stdout.write('\x1b[1A')
    sys.stdout.write('\x1b[2K')

### Variablen / Daten #########################################################

roh = 1.80 #g/ccm Massendichte vom Absorber (hier Graphit)
E1 = 0.662 #MeV Vorgegebene Emissionsenergie v. Cs-137 (eig. 661,657 keV IAEA)
re = 2.8179403262E-15 #m Elektronenradius (Wikipedia)
E0 = 0.511 #MeV Elektronenruheenergie

Teilchenanzahl = 1e5
Dicke = input('Dicke des Absorbers d ∈ N in cm: ') or 5 #muss ganzzahlig sein
Shift =input('Verschiebung des Detektors in z-Richtung (vertikal) in cm: ')or 0
N0 = int(Teilchenanzahl)
d = int(Dicke)
s = int(Shift)   #wird geändert wenn der Detektor verschoben wird in z-Richtung


l = 10-d/2  #l definiert den Anfang des Absorbers auf der x-Achse
r = 10+d/2  #r definiert das Ende des Absorbers auf der x-Achse
counter = 0 #Anzahl der Photonen emittiert in -x-Richtung
absorbiert = 0 #Anzahl der im Absorber absorbierten Photonen
durchdrung = 0 #Anzahl der Photonen, die den Absorber durchdrungen haben
reflektiert = 0 #Anzahl der Photonen, die am/im Absorber zurück gestreut werden
detektiert = 0 #Anzahl der Photonen, die die Detektorplatte treffen

Ortsauflösung = np.zeros((10,10))  #leeres Array für Ortsauflösung vom Detektor
Spektrum = np.empty(0)             #leeres Energiespektrum vom Detektor

bar = '█'
negativ_bar = '-'

df = pd.read_csv('Crossections_Carbon.dat', delim_whitespace=True, skiprows=6, 
                 header=None, names=[
    'PHOTON ENERGY (MeV)',
    'COHERENT SCATTERING (cm2/g)',
    'INCOHERENT SCATTERING (cm2/g)',
    'PHOTOELECTRIC ABSORPTION (cm2/g)',
    'PAIR PRODUCTION IN NUCLEAR FIELD (cm2/g)',
    'PAIR PRODUCTION IN ELECTRON FIELD (cm2/g)',
    'TOTAL ATTENUATION WITH COHERENT SCATTERING (cm2/g)',
    'TOTAL ATTENUATION WITHOUT COHERENT SCATTERING (cm2/g)'])

cs = np.genfromtxt('Crossections_Carbon.dat', dtype=float, skip_header=6,
                   autostrip=True)

### Berechnungen ##############################################################

for i in range(int(N0)):

    E1 = 0.662   #setzt die Energe zurrück auf die Emissionsenergie,
                 #wenn die Schleife ein neues Photon berechnet
    mu = np.random.uniform(-1,1)
    azi = np.random.uniform(0, 2*np.pi)
  
    while mu <= 0:
        mu = np.random.uniform(-1,1)
        counter += 1
        
        #die 2 if-Bedingungen brechen die for-Schleife ab, falls die 
        #Teilchenzahl erreicht wird durch Emission in negative x-Richtung
        if counter + absorbiert + reflektiert + durchdrung >= N0:
            break
    if counter + absorbiert + reflektiert + durchdrung >= N0:
        break 
   
    #Interpolation für Massenschwächungskoeffizienten von Photo-/Compton-Effekt
    cs_p = np.interp(E1, cs[:,0], cs[:,3])  #Photoeffekt
    cs_c = np.interp(E1, cs[:,0], cs[:,2])  #Compton-Effekt
    
    Sigma_p = cs_p*roh  #Berechnung Wirkungsquerschnittsdichte für Photoeffekt
    Sigma_c = cs_c*roh  #Berechnung Wirkungsquerschnittsdichte für Compton
    Sigma_tot = Sigma_p+Sigma_c

    x = np.random.uniform(0,1)     #auslosen Zufallsvariable für freie Weglänge
    y = -1/Sigma_tot * np.log(x)   #Berechnung frei Weglänge

    #Flugvektor aus Emission der Quelle, später für Streuung nach einer WW
    vektor = np.array([mu, np.sqrt(1-mu**2)*np.cos(azi), 
                       np.sqrt(1-mu**2)*np.sin(azi)])
    #Skaliert Flugvektor, sodass Photon auf die Absorber trifft
    Ortsvektor = vektor / mu * (10-d/2)
    
    #Berechnet Eindringtiefe in die Absorber mit freier Weglänge
    Ortsvektor = Ortsvektor + Ortsvektor / np.linalg.norm(Ortsvektor) * y
    
    #Abbruchbedingungen: Photone verlässt den Absorber oder 
    #gibt seine Energie vollständig ab
    while  E1 > 0.01 and Ortsvektor[0] >= l and Ortsvektor[0] <= r:
        if counter + absorbiert + reflektiert + durchdrung >= N0:
            break
    
        p_p = Sigma_p/Sigma_tot     #Wahrscheinlichkeit für Photoeffekt
        p = np.random.uniform(0,1)  #Würfelt ob Photoeffekt eintritt
        
        if p < p_p: #bricht while-Schleife ab im Fall eines Photoeffekts
              
            absorbiert += 1
            break        
            
        else:
            #Normierung des Flug-Vektors
            xf = vektor[0] / np.linalg.norm(vektor)
            yf = vektor[1] / np.linalg.norm(vektor)
            zf = vektor[2] / np.linalg.norm(vektor)
            
            #Flugvektoren in rotiertem Koordinatensystem (K-Tilde)
            if xf != 0 and yf != 0:
                j_tilde = np.array([1/xf / np.sqrt(1/xf**2+1/yf**2),
                     -1/yf / np.sqrt(1/xf**2 + 1/yf**2), 0])
            elif xf == 0 and yf != 0:
                j_tilde = np.array([1,0,0])
            elif yf == 0:
                j_tilde = np.array([0,1,0])
                
            i_tilde = np.array([xf,yf,zf])
            k_tilde = np.cross(i_tilde,j_tilde)
            
            #Rotationsmatrix zur Rücktransformation in globales System (K)
            R = np.column_stack((i_tilde, j_tilde, k_tilde))

            #auslosen von Streuwinkeln
            mu2 = Verwerfungsmethode(E1)         
            azi = np.random.uniform(0 , 2 * np.pi)
    
            #Berechnet neue Flugrichtung nach WW im rotiertem System (K-Tilde)
            x_tilde_neu = mu2
            y_tilde_neu = np.sqrt(1-mu2**2) * np.cos(azi)
            z_tilde_neu = np.sqrt(1-mu2**2) * np.sin(azi)
                              
            v_tilde = np.array([x_tilde_neu,y_tilde_neu,z_tilde_neu])
            
            #Rücktransformation neuer Flugrichtung in globales System (K)
            # vektor = np.dot(v_tilde,R) ergibt nicht das selbe wie v_tilde@R
            vektor = R @ v_tilde
            
            #Berechnung neuer Energie des Photons nach WW 
            gamma = E1 / E0
            E1 = E1 / (1+gamma*(1-mu2))
            
            #Abbruch while-Schleife falls gesamte Energie abgegeben
            if E1 < 0.01:
                absorbiert +=1
                break
            
            #Berechnet neue Massenschwächungskoeffizienten für Energie nach WW
            cs_p = np.interp(E1, cs[:,0], cs[:,3])
            cs_c = np.interp(E1, cs[:,0], cs[:,2])
            
            Sigma_p = cs_p*roh
            Sigma_c = cs_c*roh
            Sigma_tot = Sigma_p+Sigma_c
            
            #neue freie Weglänge zur nächsten WW/aus dem Absorber
            x = np.random.uniform(0,1)
            y = -1/Sigma_tot * np.log(x)   
            
            #neuer WW Ort falls im Absorber / Austritt aus dem Absorber
            #verwendet für Abbruchbedingung der while-Schleife
            Ortsvektor = Ortsvektor + vektor * y
           
    #Falls letztes Photon absorbiert wurde bricht for-Schleife ab
    if counter + absorbiert + reflektiert + durchdrung >= N0:
        break
    #Zählt Anzahl Photonen, die den Absorber durchdringen
    if Ortsvektor[0] > r:
        durchdrung += 1
        
        #verlängert den Vektor bis zu Detektorebene
        hektor = Ortsvektor + vektor * ((20-Ortsvektor[0])/vektor[0])
        
        #Überprüft y- und z-Koordinaten ob und wo diese den Detektor treffen
        if -5 <= hektor[1] < 5 and -5+s <= hektor[2] < 5+s:
            ya = int(np.floor(hektor[1])+5)
            za = int(abs(np.floor(hektor[2]+1)-5-s))
            
            Ortsauflösung[za,ya] += 1
            Spektrum = np.append(Spektrum,E1)
            detektiert += 1
            
    elif Ortsvektor[0] < l:
        reflektiert += 1
    
    progress =(absorbiert + reflektiert + durchdrung + counter) / N0
    progress_bar = bar*int(50*progress)+negativ_bar*(50-int(progress*50))
    
    sys.stdout.write(f'\r  Progress|{progress_bar}|: {progress*100:.2f}%')
    sys.stdout.flush()
    
    # print(f'Berechne : {progress} %', end='\r')
sys.stdout.write('\n')
### Auswertung / Plots ########################################################

#Auswertung der Transmissino/Reflexion/Absorption in der Konsole    
print('Anteil der den Absorber durchdringenden Photonen (p+): ',durchdrung/N0)
print('Anteil der vom Absorber reflektieren Photonen (p-): ',reflektiert/N0)
print('Anteil der vom Absorber absorbierten Photonen (p0) :',absorbiert/N0)
print('Anteil der auf den Detektor treffenden Photonen (pd): ',detektiert/N0)

#Plotted das Energiespektrum vom Detektor
plt.hist(Spektrum,bins=100,log=True)
plt.title('detektiertes Energiespektrum des Detektors')
plt.xlabel('Energie in MeV')
plt.ylabel('Anzahl Ereignisse')
plt.show()

#Plotted die Ortsauflösung vom Detektor
plt.imshow(Ortsauflösung, extent =[-5, 5, -5+s, 5+s], cmap ='viridis', vmin =0)
plt.colorbar(label ='Anzahl detektierter Photonen')
#für den Beleg wurde die vertikale Achse als y definiert, genannt hab ich diese
#im Code z, und nur für den Plot hab ich die Beschriftung "korrekt" eingefügt
plt.xlabel('z-Achse')
plt.ylabel('y-Achse')
plt.title('Ortsauflösung des Detektors')
plt.show()


