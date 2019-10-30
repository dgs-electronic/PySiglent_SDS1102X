import visa
import matplotlib
import matplotlib.pyplot as pl
import scipy.signal.windows as win
from numpy import abs, log10
from numpy.fft import fftshift, fftfreq
from scipy.fftpack import fft, next_fast_len
from matplotlib.pyplot import plot, semilogx, grid, loglog, semilogy, figure, xlabel, ylabel, legend, show



matplotlib.use('TkAgg')


# Modul zur Steuerung eines Siglent SDS1102X DSO
class Trace:
    def __init__(self):
        self.time = []
        self.volt = []
        self.freq = []
        self.spectrum = []
        self.samplingrate = 0
        self.dso = 0
        self.timedivision = 0
        self.voltagedivision = 0
        self.offset = 0

    def getIdentity(self):
        return self.dso.query("*IDN?")

    def getAcquisitionMode(self):
        return self.dso.query("*ACQW?")

    def getProbeAttenuation(self, channel):
        if channel == 1:
            return self.dso.query("C1:ATTN?")
        elif channel == 2:
            return self.dso.query("C2:ATTN?")

    def getVoltageDivision(self,channel):
        return self.dso.query("c" + str(channel) + ":vdiv?")

    def getVoltageOffset(self,channel):
        return self.dso.query("c" + str(channel) + ":ofst?")

    def getTimeDivision(self):
        return self.dso.query("tdiv?")

    def getSamplingRate(self):
        sara = self.dso.query("sara?")
        sara_unit = {'G': 1E9, 'M': 1E6, 'k': 1E3, 'K': 1E3}
        for unit in sara_unit.keys():
            if sara.find(unit) != -1:
                sara = sara.split(unit)
                sara = float(sara[0]) * sara_unit[unit]
                break
        sara = float(sara)
        return sara

    def getTrace(self,channel):
        self.dso.write("chdr off")                #Keine Befehlsheader nur die Antwort zeigen
        self.dso.timeout = 60000                  #60 Sekunden Timeout
        self.dso.chunk_size = 20 * 1024 * 1024
        vdiv = self.getVoltageDivision(channel)   #Spannungseinheit holen
        ofst = self.getVoltageOffset(channel)     #Spannungsoffset holen
        tdiv = self.getTimeDivision()             #Zeitbasis holen
        sara = self.getSamplingRate()             #Abtastrate holen
        self.dso.write("c" + str(channel) + ":wf? dat2")
        recv = list(self.dso.read_raw())[15:]
        recv.pop()
        recv.pop()
        #Spannungs Rohwerte umrechnen
        volt_value = []
        for data in recv:
            if data > 127:
                data = data - 255
            else:
                pass
            volt_value.append(data)
        #Zeit Rohwerte umrechnen
        time_value = []
        for idx in range(0, len(volt_value)):
            volt_value[idx] = volt_value[idx] / 25 * float(vdiv) - float(ofst)
            time_data = -(float(tdiv) * 14 / 2) + idx * (1 / sara)
            time_value.append(time_data)
        #Werte in Objekt-Variablen überspielen
        self.voltagedivision = vdiv
        self.timedivision = tdiv
        self.samplingrate = sara
        self.offset = ofst
        self.time = time_value
        self.volt = volt_value
        return '-DONE-'

    def calcFFT(self):
        #Länge des zu transformierenden Vectors bestimmen
        N = len(self.volt)
        #Optimale FFT-Länge bestimmen und ggf Zeropadding
        N_fast = next_fast_len(N)
        #Kaiser-Bessel-Fenster bestimmen
        w = win.kaiser(N,12)
        w = w / w.sum() * N #Normierung damit Amplituden wieder stimmen
        #Transformieren
        yf = fft(self.volt*w, N_fast)
        yf = fftshift(yf)
        #Abschneiden
        yf = yf[N_fast//2:]
        #Skalieren
        yf = 2/N_fast*abs(yf)
        #Frequenzachse erzeugen
        xf = fftfreq(N_fast)
        xf = fftshift(xf)
        #Abschneiden
        xf = xf[N_fast//2:]
        #Skalieren
        xf = xf*self.samplingrate
        #Veröffentlichen der Werte
        self.spectrum = yf
        self.freq = xf
        return '-DONE-'

    def plotSpectrum(self,titel="Spectrum"):
        figure(figsize=(7, 5))
        semilogx(self.freq, self.spectrum, 'r', markersize=1, label=titel)
        xlabel("Frequency [Hz]")
        ylabel("Amplitude [V]")
        legend()
        grid(b=True, which='major', color='g', linestyle='-')
        grid(b=True, which='minor', color='g', linestyle='-')
        show()
        return '-DONE-'

    def plotPowerSpectrum_dBm(self,titel="PowerSpectrum", R_bezug=1):
        spectrum_dB = 10*log10(self.spectrum*self.spectrum / 1E-3 / R_bezug)
        figure(figsize=(7, 5))
        semilogx(self.freq, spectrum_dB, 'r', markersize=1, label=titel)
        xlabel("Frequency [Hz]")
        ylabel("Amplitude [dBm]")
        legend()
        grid(b=True, which='major', color='g', linestyle='-')
        grid(b=True, which='minor', color='g', linestyle='-')
        show()
        return '-DONE-'

    def plotPowerSpectrum(self,titel="PowerSpectrum"):
        spectrum_pw = (self.spectrum * self.spectrum)
        figure(figsize=(7, 5))
        semilogx(self.freq, spectrum_pw, 'r', markersize=1, label=titel)
        xlabel("Frequency [Hz]")
        ylabel("Amplitude [V²]")
        legend()
        grid(b=True, which='major', color='g', linestyle='-')
        grid(b=True, which='minor', color='g', linestyle='-')
        show()
        return '-DONE-'

    def plotTrace(self,titel="Trace1"):
        figure(figsize=(7, 5))
        plot(self.time, self.volt, 'r', markersize=1, label=titel)
        xlabel("Frequency [s]")
        ylabel("Amplitude [V]")
        legend()
        grid(b=True, which='major', color='g', linestyle='-')
        grid(b=True, which='minor', color='g', linestyle='-')
        show()
        return '-DONE-'







# main
_rm = visa.ResourceManager()
sds = _rm.open_resource("TCPIP::192.168.178.66::INSTR")

TestScope = Trace()
TestScope.dso = sds
print(TestScope.getIdentity())
print(TestScope.getTrace(1))
print(TestScope.calcFFT())
print(TestScope.plotSpectrum("Testspektrum"))
#FFT-Berechnung
#N = len(TestScope.volt)
#N_fast = next_fast_len(N)
#w = win.kaiser(N,12)
#w = w / w.sum() * N

#test = fft(TestScope.volt*w,N_fast)
#test = fftshift(test)
#test = test[N_fast//2:]
#test = 2/N_fast*abs(test)
#test_dBv = 20 * log10(test)
#test_dBm = 10 * log10(test * test / 50 / 1E-3)

#xf = fftfreq(N_fast)
#xf = fftshift(xf)
#xf = xf[N_fast//2:]
#xf = xf*TestScope.samplingrate
#FFT-Berechnung ENDE


#pl.figure(figsize=(7, 5))
#pl.plot(TestScope.time, TestScope.volt, 'g', markersize=2, label=u"Spannung an C337")
#pl.semilogx(TestScope.freq, TestScope.spectrum, 'r', markersize=1, label=u"Spannung an C337")
#pl.xlabel("Time in [µs]")
#pl.ylabel("Voltage in [V]")
#pl.legend()
#pl.grid()
#pl.grid(b=True, which='major', color='g', linestyle='-')
#pl.grid(b=True, which='minor', color='g', linestyle='-')
#pl.show()





