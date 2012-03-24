#! /usr/bin/env python

import os
import sys
import math
import string


def box(pos, hsz=[1,1,1], material='light gray', ori=None, vel=None, angvel=None):
  s = '<box hx="%f" hy="%f" hz="%f">\n' % (hsz[0], hsz[1], hsz[2]);
  s += '  <pos x="%f" y="%f" z="%f"/>\n' % (pos[0], pos[1], pos[2]);
  s += '  <bodymaterial name="%s" />\n' % material;
  if (ori != None):
    s += '  <ori x="%f" y="%f" z="%f" theta="%f"/>\n' % (ori[0], ori[1], ori[2], ori[3]);
  if (vel != None):
    s += '  <vel x="%f" y="%f" z="%f"/>\n' % (vel[0], vel[1], vel[2]);
  if (angvel != None):
    s += '  <angvel x="%f" y="%f" z="%f"/>\n' % (angvel[0], angvel[1], angvel[2]);
  s += '</box>';
  return s;

def sphere(pos, r=1, material='dark gray', ori=None, vel=None, angvel=None):
  s = '<sphere r="%f">\n' % (r);
  s += '  <pos x="%f" y="%f" z="%f"/>\n' % (pos[0], pos[1], pos[2]);
  s += '  <bodymaterial name="%s" />\n' % material;
  if (ori != None):
    s += '  <ori x="%f" y="%f" z="%f" theta="%f"/>\n' % (ori[0], ori[1], ori[2], ori[3]);
  if (vel != None):
    s += '  <vel x="%f" y="%f" z="%f"/>\n' % (vel[0], vel[1], vel[2]);
  if (angvel != None):
    s += '  <angvel x="%f" y="%f" z="%f"/>\n' % (angvel[0], angvel[1], angvel[2]);
  s += '</sphere>';
  return s;


if (__name__ == '__main__'):
  print(box([1,2,3]));
  print(sphere([1,2,3]));

