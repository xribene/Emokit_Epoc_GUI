# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 04:50:40 2017

@author: Christodoulos Benetatos - xribene
"""
###############################################################################
from __future__ import division
from pyqtgraph.Qt import QtCore, QtGui
from emokit.emotiv import Emotiv
import matplotlib.pyplot as plt
from collections import deque
from scipy import signal,fft
from Queue import Queue
import time, threading
import pyqtgraph as pg
import numpy as np
import sys
###############################################################################
# Helper functions

def eeg_fft(y,Fs=128,show=False,limits=[0,30,0,20]) : 
    y=np.atleast_2d(y)
    [C,N]=np.shape(y)
    T=1/Fs
    y_f = fft(y)
    xf = np.linspace(0.0, 1.0/(2.0*T), N/2,endpoint=False)
    xf=np.atleast_2d(xf)
    xf=np.tile(xf,(C,1))
    yf= (np.square(2.0/N *np.abs(y_f[:,0:N//2])))
    if show:
        plt.ion()
        plt.plot(np.squeeze(xf), np.squeeze(yf))
        plt.axis(limits)
        plt.xticks(range(limits[1]))
        plt.grid()
    return yf,xf

def filtering(eeg,cut_off=[2,40],mode='band',order=4,show=False,limits=[0,30,0,20],ex=[111,111]):
    eeg=np.atleast_2d(eeg)
    [C,N]=np.shape(eeg)
    tmp=np.zeros(shape=(np.shape(eeg)))
    b, a = signal.butter(order, np.array((cut_off))/(Fs/2), btype = mode)
    for i in xrange(C):
        if (i==ex[0]) or (i==ex[1]):
            tmp[i,:]=eeg[i,:]
            continue
        tmp[i,:]=signal.filtfilt(b, a, eeg[i,:] ) 
        #print np.mean(tmp[i,:])
    if show:
        eeg_fft(np.squeeze(eeg),Fs,show,limits)
        eeg_fft(np.squeeze(tmp),Fs,show,limits)
    return tmp

def quality_color(av):
    aa=av//20
    aa=(aa<255)*aa+(aa>255)*255+(aa==255)*255
    return (255-aa, aa, 0)
def pad(array_in, result):
    # zero pad array so that shape(array_in)=shape(result)
    [a,b]=np.shape(array_in)
    [k,l]=np.shape(result)
    start1=int((k-a)/2)
    start2=int((l-b)/2)
    result[start1:(start1+a),start2:(start2+b)]=array_in
    return result
def next_pow(x):
    return 1<<(x-1).bit_length()
##############################################################################
class ring_buffer(object):
    def __init__(self,size):
        self.size=size
        self._buffered= deque([], self.size)
    def write(self, value):
        self._buffered.append(value)
    def write_ex(self, value):
        self._buffered.extend(value)
    def show(self):
        print(self._buffered)
    def copy(self,overlap):
        tmp=list(self._buffered)
        return tmp[0:overlap]
    def calls(self):
        return self.write.calls
    def list_ret(self):
        a=list(self._buffered)
        return a

class Plotter():
    def __init__(self,electrodes,tw_sec,step,q1,flag1):
        self.q1=q1
        self.flag1=flag1
        self.step=step
        self.channels=len(electrodes.split(' '))
        self.curve=[]
        self.p=[]
        self.scores=[]
        self.now=0
        self.now2=0
        self.ptr1=0
        self.Fs=128
        self.N=tw_sec*self.Fs
        # array with reference shape for zeropad in line 146
        self.b=np.zeros(shape=(self.channels,next_pow(self.N))) 
        if show:
            self.win = pg.GraphicsWindow()
            self.win.setWindowTitle('Emotiv Epoc EEG Data')
            for i in xrange(2*self.channels):
                if i<self.channels:
                    self.p.append(self.win.addPlot(colspan=self.channels,title=electrodes.split(' ')[i]))
                    self.win.nextRow()
                if i>(self.channels-1):
                    self.p.append(self.win.addPlot(colspan=1,title=electrodes.split(' ')[i-self.channels]+' PSD'))
            
            data1 = np.random.normal(size=10)
            for i in xrange(2*self.channels):
                self.curve.append(self.p[i].plot(data1))
            self.text_peak=[]
            self.arrow=[]
            self.text_qual=[]
            for i in xrange(2*self.channels):
                if i<self.channels:
                    #p[i].setXRange(0,1000, padding=0)
                    self.p[i].setYRange(-200, 200, padding=0)
                    tmp_text=pg.TextItem(anchor=(-0.4,1.6), fill=(0, 0, 255, 100))
                    self.text_qual.append(tmp_text)
                    self.p[i].addItem(tmp_text)
                if i>(self.channels-1):
                    self.p[i].setXRange(0, 50, padding=0)
                    self.p[i].setYRange(0, 20, padding=0)
                    tmp_text=pg.TextItem(anchor=(0.5,2), fill=(0,0,255, 80))
                    tmp_arrow=pg.ArrowItem( angle=-90,brush=(0,0,255))
                    self.p[i].addItem(tmp_text)
                    self.p[i].addItem(tmp_arrow)
                    self.text_peak.append(tmp_text)
                    self.arrow.append(tmp_arrow)
        self.set_timer()
        
    def update_plots(self):
        #print time.time()-self.now2
        #self.now2=time.time()
        y=np.array(self.q1.get(block=True, timeout=5)).transpose(1,0,2)
        qual=y[:,:,1]
        val_a=y[:,:,0]
        val=filtering(val_a,cut_off=[2 ,30],mode='band')
        val_pad=pad(val,self.b)
        valf,xf=eeg_fft(val_pad,self.Fs,show=False,limits=[0,30,0,20])
        self.ptr1 += self.step
        if show:
            for i in xrange(2*self.channels):
                if i<self.channels:
                    av_qual=qual[i].mean()
                    qual_color=quality_color(av_qual)
                    self.curve[i].setData(val[i],pen=pg.mkPen(color=qual_color),width=2)
                    self.curve[i].setPos(self.ptr1, 0)
                    self.text_qual[i].setText('%0.1f' % av_qual)
                    self.text_qual[i].setColor(color=(0,255,255))
                    self.text_qual[i].setPos(self.ptr1,0)
                if i>(self.channels-1):
                    self.curve[i].setData(xf[0],(valf[i-self.channels])) 
                    self.text_peak[i-self.channels].setText('%0.3f' % xf[0][np.argmax(valf[i-self.channels])])
                    self.text_peak[i-self.channels].setColor(color=(0,255,255))
                    self.text_peak[i-self.channels].setPos(xf[0][np.argmax(valf[i-self.channels])], valf[i-self.channels].max())
                    self.arrow[i-self.channels].setPos(xf[0][np.argmax(valf[i-self.channels])], valf[i-self.channels].max())
    def check_flag(self):
        self.flag1.wait()
        #self.now=time.time()
        self.update_plots()
        #print time.time()-self.now
        pass
    def set_timer(self):
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.check_flag)
        self.timer.start(10)
    
# reads data from emotiv and sends them to Plotter every 'step/Fs' seconds,
# through q1
class Reader(threading.Thread):
    def __init__(self,q1,flag1,step,tw,electrodes):
        super(Reader, self).__init__()
        self.Fs=128
        self._stop = threading.Event()
        self.q1=q1
        self.flag1=flag1
        self.step=step
        self.tw=tw*self.Fs
        self.electrodes=electrodes
    def stop(self):
        self._stop.set()
    def stopped(self):
        return self._stop.isSet()
    def run(self):
        O1_buff = ring_buffer(self.tw)
        i=0
        #old=0
        #old2=0
        with Emotiv(display_output=False, verbose=True, write=False) as headset:
            try:
                while not self._stop.isSet():
                    packet = headset.dequeue()
                    if packet is not None:
                        i=i+1
                        data=[]
                        for name in electrodes.split(' '):
                            data.append([packet.sensors[name]['value'],packet.sensors[name]['quality']])
                        O1_buff.write(data) 
                        #print(time.time()-old)
                        #old=time.time()
                        pass
                    if i==self.step:
                        self.flag1.set()
                        self.q1.put(O1_buff.list_ret())
                        self.flag1.clear()
                        i=0
                        #print(time.time()-old2)
                        #old2=time.time()
            except :
                pass

##############################################################################
if __name__ == '__main__':
    Fs=128
    # O1 O2 P7 P8 AF3 F7 F3 FC5 T7 T8 FC6 F4 F8 AF4 X Y
    electrodes='O1 O2' # choose which sensors to graph
    tw_sec=2 # time window in which fft will be calculated 
    q1=Queue()
    step=np.round(0.5*Fs) # how many new points will be graphed in every update
    # or how many seconds (0.5) between update_plots repetetions
    # or overlap between tw_sec windows
    flag1=threading.Event()
    thread1=Reader(q1,flag1,step,tw_sec,electrodes)
    show=1
    app = QtGui.QApplication(sys.argv)
    s = Plotter(electrodes,tw_sec,step,q1,flag1)
    thread1.start()
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_() 
    thread1.stop()
