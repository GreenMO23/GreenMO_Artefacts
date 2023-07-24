#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: X300 Fdma
# GNU Radio version: 3.7.14.0
##################################################

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

from PyQt4 import Qt
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import pmt
import sys
import time
from gnuradio import qtgui


class x300_fdma(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "X300 Fdma")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("X300 Fdma")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "x300_fdma")
        self.restoreGeometry(self.settings.value("geometry").toByteArray())


        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 10e6
        self.gain2 = gain2 = 25
        self.gain1 = gain1 = 20
        self.center_freq_3 = center_freq_3 = 2.462e9
        self.center_freq_2 = center_freq_2 = 2.442e9
        self.center_freq_1 = center_freq_1 = 2.452e9
        self.center_freq_0 = center_freq_0 = 2.432e9

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_sink_0 = uhd.usrp_sink(
        	",".join(("addr0=192.168.10.2,addr1=192.168.11.2", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(4),
        	),
        )
        self.uhd_usrp_sink_0.set_clock_source('external', 0)
        self.uhd_usrp_sink_0.set_time_source('external', 0)
        self.uhd_usrp_sink_0.set_subdev_spec('A:0 B:0', 0)
        self.uhd_usrp_sink_0.set_clock_source('external', 1)
        self.uhd_usrp_sink_0.set_time_source('external', 1)
        self.uhd_usrp_sink_0.set_subdev_spec('A:0 B:0', 1)
        self.uhd_usrp_sink_0.set_samp_rate(samp_rate)
        self.uhd_usrp_sink_0.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_sink_0.set_center_freq(center_freq_0, 0)
        self.uhd_usrp_sink_0.set_gain(gain2, 0)
        self.uhd_usrp_sink_0.set_antenna('TX/RX', 0)
        self.uhd_usrp_sink_0.set_bandwidth(samp_rate, 0)
        self.uhd_usrp_sink_0.set_center_freq(center_freq_1, 1)
        self.uhd_usrp_sink_0.set_gain(gain2, 1)
        self.uhd_usrp_sink_0.set_antenna('TX/RX', 1)
        self.uhd_usrp_sink_0.set_bandwidth(samp_rate, 1)
        self.uhd_usrp_sink_0.set_center_freq(center_freq_2, 2)
        self.uhd_usrp_sink_0.set_gain(gain2, 2)
        self.uhd_usrp_sink_0.set_antenna('TX/RX', 2)
        self.uhd_usrp_sink_0.set_bandwidth(samp_rate, 2)
        self.uhd_usrp_sink_0.set_center_freq(center_freq_3, 3)
        self.uhd_usrp_sink_0.set_gain(gain2, 3)
        self.uhd_usrp_sink_0.set_antenna('TX/RX', 3)
        self.uhd_usrp_sink_0.set_bandwidth(samp_rate, 3)
        self.blocks_file_source_0 = blocks.file_source(gr.sizeof_gr_complex*1, '/mnt/ramdisk/uhd_ofdm_tx_1user_coded_QAM16_34r_0.dat', True)
        self.blocks_file_source_0.set_begin_tag(pmt.PMT_NIL)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_file_source_0, 0), (self.uhd_usrp_sink_0, 0))
        self.connect((self.blocks_file_source_0, 0), (self.uhd_usrp_sink_0, 1))
        self.connect((self.blocks_file_source_0, 0), (self.uhd_usrp_sink_0, 2))
        self.connect((self.blocks_file_source_0, 0), (self.uhd_usrp_sink_0, 3))

    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "x300_fdma")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_sink_0.set_samp_rate(self.samp_rate)
        self.uhd_usrp_sink_0.set_bandwidth(self.samp_rate, 0)
        self.uhd_usrp_sink_0.set_bandwidth(self.samp_rate, 1)
        self.uhd_usrp_sink_0.set_bandwidth(self.samp_rate, 2)
        self.uhd_usrp_sink_0.set_bandwidth(self.samp_rate, 3)

    def get_gain2(self):
        return self.gain2

    def set_gain2(self, gain2):
        self.gain2 = gain2
        self.uhd_usrp_sink_0.set_gain(self.gain2, 0)

        self.uhd_usrp_sink_0.set_gain(self.gain2, 1)

        self.uhd_usrp_sink_0.set_gain(self.gain2, 2)

        self.uhd_usrp_sink_0.set_gain(self.gain2, 3)


    def get_gain1(self):
        return self.gain1

    def set_gain1(self, gain1):
        self.gain1 = gain1

    def get_center_freq_3(self):
        return self.center_freq_3

    def set_center_freq_3(self, center_freq_3):
        self.center_freq_3 = center_freq_3
        self.uhd_usrp_sink_0.set_center_freq(self.center_freq_3, 3)

    def get_center_freq_2(self):
        return self.center_freq_2

    def set_center_freq_2(self, center_freq_2):
        self.center_freq_2 = center_freq_2
        self.uhd_usrp_sink_0.set_center_freq(self.center_freq_2, 2)

    def get_center_freq_1(self):
        return self.center_freq_1

    def set_center_freq_1(self, center_freq_1):
        self.center_freq_1 = center_freq_1
        self.uhd_usrp_sink_0.set_center_freq(self.center_freq_1, 1)

    def get_center_freq_0(self):
        return self.center_freq_0

    def set_center_freq_0(self, center_freq_0):
        self.center_freq_0 = center_freq_0
        self.uhd_usrp_sink_0.set_center_freq(self.center_freq_0, 0)


def main(top_block_cls=x300_fdma, options=None):

    from distutils.version import StrictVersion
    if StrictVersion(Qt.qVersion()) >= StrictVersion("4.5.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()
    tb.start()
    tb.show()

    def quitting():
        tb.stop()
        tb.wait()
    qapp.connect(qapp, Qt.SIGNAL("aboutToQuit()"), quitting)
    qapp.exec_()


if __name__ == '__main__':
    main()
