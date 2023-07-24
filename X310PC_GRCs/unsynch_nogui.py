#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Unsynch Nogui
# GNU Radio version: 3.7.14.0
##################################################

from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import pmt
import time


class unsynch_nogui(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Unsynch Nogui")

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 10e6
        self.head_samps = head_samps = 1200e6
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
        self.uhd_usrp_sink_0.set_time_now(uhd.time_spec(time.time()), uhd.ALL_MBOARDS)
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
        self.blocks_head_0_0_0_0 = blocks.head(gr.sizeof_gr_complex*1, int(head_samps))
        self.blocks_head_0_0_0 = blocks.head(gr.sizeof_gr_complex*1, int(head_samps))
        self.blocks_head_0_0 = blocks.head(gr.sizeof_gr_complex*1, int(head_samps))
        self.blocks_head_0 = blocks.head(gr.sizeof_gr_complex*1, int(head_samps))
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
        self.connect((self.blocks_file_source_0_0_0_0, 0), (self.blocks_head_0_0_0_0, 0))
        self.connect((self.blocks_file_source_0_2, 0), (self.blocks_head_0, 0))
        self.connect((self.blocks_file_source_0_2_0, 0), (self.blocks_head_0_0, 0))
        self.connect((self.blocks_file_source_0_2_0_0, 0), (self.blocks_head_0_0_0, 0))
        self.connect((self.blocks_head_0, 0), (self.uhd_usrp_sink_0, 0))
        self.connect((self.blocks_head_0_0, 0), (self.uhd_usrp_sink_0, 1))
        self.connect((self.blocks_head_0_0_0, 0), (self.uhd_usrp_sink_0, 2))
        self.connect((self.blocks_head_0_0_0_0, 0), (self.uhd_usrp_sink_0, 3))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_sink_0.set_samp_rate(self.samp_rate)
        self.uhd_usrp_sink_0.set_bandwidth(self.samp_rate, 0)
        self.uhd_usrp_sink_0.set_bandwidth(self.samp_rate, 1)
        self.uhd_usrp_sink_0.set_bandwidth(self.samp_rate, 2)
        self.uhd_usrp_sink_0.set_bandwidth(self.samp_rate, 3)

    def get_head_samps(self):
        return self.head_samps

    def set_head_samps(self, head_samps):
        self.head_samps = head_samps
        self.blocks_head_0_0_0_0.set_length(int(self.head_samps))
        self.blocks_head_0_0_0.set_length(int(self.head_samps))
        self.blocks_head_0_0.set_length(int(self.head_samps))
        self.blocks_head_0.set_length(int(self.head_samps))

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


def main(top_block_cls=unsynch_nogui, options=None):

    tb = top_block_cls()
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()
