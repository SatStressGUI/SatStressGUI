#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 09:12:26 2017

@author: ndoan
"""
import shutil 
class A: 
    variable1 = True 
    
    def __init__(self): 
        print("hello, an instance of class a was created")

class B: 
    def __init__(self): 
        print("an instance of class b was created")
        print(A.variable1)
        
b = B() 