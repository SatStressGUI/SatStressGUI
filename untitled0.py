#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 16:31:45 2017

@author: ndoan
"""

import sys
import datetime
import imageio

VALID_EXTENSIONS = ('png', 'jpg')


def create_gif(filenames, duration):
    images = []
    for filename in filenames:
        im = imageio.imread(filename) 
        #im.shape(200,200,2)
        images.append(im)
    output_file = 'Gif-%s.gif' % datetime.datetime.now().strftime('%Y-%M-%d-%H-%M-%S')
    imageio.mimsave(output_file, images, duration=duration)


if __name__ == "__main__":
    duration = .5
    filenames = ["/Users/ndoan/Desktop/orbit_360.00.png", "/Users/ndoan/Desktop/orbit_324.00.png","/Users/ndoan/Desktop/orbit_288.00.png"]

    create_gif(filenames, duration)