#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: X300 Twinrx
# GNU Radio version: 3.7.14.0
##################################################

from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import time
import os

class x300_twinrx(gr.top_block):

    def __init__(self,test_id,test_str):
        gr.top_block.__init__(self, "X300 Twinrx")

        ##################################################
        # Variables
        ##################################################
        self.test_str = test_str 
        self.test_id = test_id 
        self.samp_rate_tx = samp_rate_tx = 12.5e6
        self.samp_rate = samp_rate = 50e6
        self.rx_gain = rx_gain = 70
        self.count = count = int(2e6)
        self.center_freq = center_freq = 2.484e9

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_source_0 = uhd.usrp_source(
        	",".join(("addr=192.168.40.2", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(2),
        	),
        )


        

        self.uhd_usrp_source_0.set_subdev_spec('A:0 A:1', 0)
        self.uhd_usrp_source_0.set_samp_rate(samp_rate)
        self.uhd_usrp_source_0.set_center_freq(center_freq, 0)
        self.uhd_usrp_source_0.set_gain(rx_gain, 0)
        self.uhd_usrp_source_0.set_antenna('RX1', 0)
        self.uhd_usrp_source_0.set_bandwidth(samp_rate, 0)
        self.uhd_usrp_source_0.set_auto_dc_offset(True, 0)
        self.uhd_usrp_source_0.set_auto_iq_balance(True, 0)
        self.uhd_usrp_source_0.set_center_freq(center_freq, 1)
        self.uhd_usrp_source_0.set_gain(rx_gain, 1)
        self.uhd_usrp_source_0.set_antenna('RX2', 1)
        self.uhd_usrp_source_0.set_bandwidth(samp_rate, 1)
        self.uhd_usrp_source_0.set_auto_dc_offset(True, 1)
        self.uhd_usrp_source_0.set_auto_iq_balance(True, 1)
        
        # Synchronizing all clocks to next pps
        self.uhd_usrp_source_0.set_time_next_pps(uhd.time_spec(0))
        time.sleep(1.1)
        self.uhd_usrp_source_0.set_start_time(uhd.time_spec_t(1.2))
        # Locking PLLs together
        # self.uhd_usrp_source_0.set_command_time(uhd.time_spec_t(1),0)
        # self.uhd_usrp_source_0.clear_command_time(0)
        # self.uhd_usrp_source_0.set_command_time(uhd.time_spec_t(1),1)
        # self.uhd_usrp_source_0.set_center_freq(center_freq, 1)
        # self.uhd_usrp_source_0.clear_command_time(1)

        self.blocks_head_0_0 = blocks.head(gr.sizeof_gr_complex*1, count)
        self.blocks_head_0 = blocks.head(gr.sizeof_gr_complex*1, count)
        self.blocks_file_sink_0_0 = blocks.file_sink(gr.sizeof_gr_complex*1, "/mnt/ramdisk/rx1_"+test_str+str(test_id)+".dat", False)
        self.blocks_file_sink_0_0.set_unbuffered(False)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_gr_complex*1, "/mnt/ramdisk/rx0_"+test_str+str(test_id)+".dat", False)
        self.blocks_file_sink_0.set_unbuffered(False)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_head_0, 0), (self.blocks_file_sink_0, 0))
        self.connect((self.blocks_head_0_0, 0), (self.blocks_file_sink_0_0, 0))
        self.connect((self.uhd_usrp_source_0, 0), (self.blocks_head_0, 0))
        self.connect((self.uhd_usrp_source_0, 1), (self.blocks_head_0_0, 0))

    def get_test_str(self):
        return self.test_str

    def set_test_str(self, test_str):
        self.test_str = test_str
        self.blocks_file_sink_0_0.open("/mnt/ramdisk/rx1_"+self.test_str+str(self.test_id)+".dat")
        self.blocks_file_sink_0.open("/mnt/ramdisk/rx0_"+self.test_str+str(self.test_id)+".dat")

    def get_test_id(self):
        return self.test_id

    def set_test_id(self, test_id):
        self.test_id = test_id
        self.blocks_file_sink_0_0.open("/mnt/ramdisk/rx1_"+self.test_str+str(self.test_id)+".dat")
        self.blocks_file_sink_0.open("/mnt/ramdisk/rx0_"+self.test_str+str(self.test_id)+".dat")

    def get_samp_rate_tx(self):
        return self.samp_rate_tx

    def set_samp_rate_tx(self, samp_rate_tx):
        self.samp_rate_tx = samp_rate_tx
        self.uhd_usrp_source_0.set_bandwidth(self.samp_rate_tx, 2)

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)
        self.uhd_usrp_source_0.set_bandwidth(self.samp_rate, 0)
        self.uhd_usrp_source_0.set_bandwidth(self.samp_rate, 1)

    def get_rx_gain(self):
        return self.rx_gain

    def set_rx_gain(self, rx_gain):
        self.rx_gain = rx_gain
        self.uhd_usrp_source_0.set_gain(self.rx_gain, 0)

        self.uhd_usrp_source_0.set_gain(self.rx_gain, 1)

        self.uhd_usrp_source_0.set_gain(self.rx_gain, 2)


    def get_count(self):
        return self.count

    def set_count(self, count):
        self.count = count
        self.blocks_head_0_0.set_length(self.count)
        self.blocks_head_0.set_length(self.count)

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.uhd_usrp_source_0.set_center_freq(self.center_freq, 0)
        self.uhd_usrp_source_0.set_center_freq(self.center_freq, 1)
        self.uhd_usrp_source_0.set_center_freq(self.center_freq, 2)


def main(top_block_cls=x300_twinrx, options=None):

    date_str = "07_20"
    folder_id = 0
    num_tests = 1
    test_str = "repeat_check"
    # sav_dir = "/home/wcsng-20/ag_data/repeat_check_"+date_str+"_"+str(folder_id)+"/"
    # os.system("mkdir "+sav_dir)

    for i in range(num_tests):
        tb = top_block_cls(i,test_str)
        tb.start()
        # try:
        #     raw_input('Press Enter to quit: ')
        # except EOFError:
        #     pass
        tb.stop()
        tb.wait()
        time.sleep(1)
        del tb
        # os.system("mv /mnt/ramdisk/rx0_"+test_str+str(i)+".dat "+sav_dir)
        # os.system("mv /mnt/ramdisk/rx1_"+test_str+str(i)+".dat "+sav_dir)
        print("Finished test #"+str(i)+", Starting next test..")
        time.sleep(1)

if __name__ == '__main__':
    main()
