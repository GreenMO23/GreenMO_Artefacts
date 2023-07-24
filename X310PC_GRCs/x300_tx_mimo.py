#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: X300 Tx Mimo
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


class x300_tx_mimo(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "X300 Tx Mimo")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("X300 Tx Mimo")
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

        self.settings = Qt.QSettings("GNU Radio", "x300_tx_mimo")
        self.restoreGeometry(self.settings.value("geometry").toByteArray())


        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 10e6
        self.gain2 = gain2 = 25
        self.gain1 = gain1 = 20
        self.file_base_str = file_base_str = "/mnt/ramdisk/uhd_ofdm_tx_4user_coded_QAM16_12r"
        self.center_freq = center_freq = 2.432e9

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
        self.uhd_usrp_sink_0.set_time_source('external', 0)
        self.uhd_usrp_sink_0.set_subdev_spec('A:0 B:0', 0)
        self.uhd_usrp_sink_0.set_time_source('external', 1)
        self.uhd_usrp_sink_0.set_subdev_spec('A:0 B:0', 1)
        self.uhd_usrp_sink_0.set_samp_rate(samp_rate)
        self.uhd_usrp_sink_0.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_sink_0.set_center_freq(center_freq, 0)
        self.uhd_usrp_sink_0.set_gain(gain2, 0)
        self.uhd_usrp_sink_0.set_antenna('TX/RX', 0)
        self.uhd_usrp_sink_0.set_bandwidth(samp_rate, 0)
        self.uhd_usrp_sink_0.set_center_freq(center_freq, 1)
        self.uhd_usrp_sink_0.set_gain(gain2, 1)
        self.uhd_usrp_sink_0.set_antenna('TX/RX', 1)
        self.uhd_usrp_sink_0.set_bandwidth(samp_rate, 1)
        self.uhd_usrp_sink_0.set_center_freq(center_freq, 2)
        self.uhd_usrp_sink_0.set_gain(gain2, 2)
        self.uhd_usrp_sink_0.set_antenna('TX/RX', 2)
        self.uhd_usrp_sink_0.set_bandwidth(samp_rate, 2)
        self.uhd_usrp_sink_0.set_center_freq(center_freq, 3)
        self.uhd_usrp_sink_0.set_gain(gain2, 3)
        self.uhd_usrp_sink_0.set_antenna('TX/RX', 3)
        self.uhd_usrp_sink_0.set_bandwidth(samp_rate, 3)
        self.blocks_file_source_0_2_0_0 = blocks.file_source(gr.sizeof_gr_complex*1, file_base_str+"_2.dat", True)
        self.blocks_file_source_0_2_0_0.set_begin_tag(pmt.PMT_NIL)
        self.blocks_file_source_0_2_0 = blocks.file_source(gr.sizeof_gr_complex*1, file_base_str+"_1.dat", True)
        self.blocks_file_source_0_2_0.set_begin_tag(pmt.PMT_NIL)
        self.blocks_file_source_0_2 = blocks.file_source(gr.sizeof_gr_complex*1, file_base_str+"_0.dat", True)
        self.blocks_file_source_0_2.set_begin_tag(pmt.PMT_NIL)
        self.blocks_file_source_0_0_0_0 = blocks.file_source(gr.sizeof_gr_complex*1, file_base_str+"_3.dat", True)
        self.blocks_file_source_0_0_0_0.set_begin_tag(pmt.PMT_NIL)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_file_source_0_0_0_0, 0), (self.uhd_usrp_sink_0, 3))
        self.connect((self.blocks_file_source_0_2, 0), (self.uhd_usrp_sink_0, 0))
        self.connect((self.blocks_file_source_0_2_0, 0), (self.uhd_usrp_sink_0, 1))
        self.connect((self.blocks_file_source_0_2_0_0, 0), (self.uhd_usrp_sink_0, 2))

    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "x300_tx_mimo")
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

    def get_file_base_str(self):
        return self.file_base_str

    def set_file_base_str(self, file_base_str):
        self.file_base_str = file_base_str
        self.blocks_file_source_0_2_0_0.open(self.file_base_str+"_2.dat", True)
        self.blocks_file_source_0_2_0.open(self.file_base_str+"_1.dat", True)
        self.blocks_file_source_0_2.open(self.file_base_str+"_0.dat", True)
        self.blocks_file_source_0_0_0_0.open(self.file_base_str+"_3.dat", True)

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.uhd_usrp_sink_0.set_center_freq(self.center_freq, 0)
        self.uhd_usrp_sink_0.set_center_freq(self.center_freq, 1)
        self.uhd_usrp_sink_0.set_center_freq(self.center_freq, 2)
        self.uhd_usrp_sink_0.set_center_freq(self.center_freq, 3)


def main(top_block_cls=x300_tx_mimo, options=None):

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
