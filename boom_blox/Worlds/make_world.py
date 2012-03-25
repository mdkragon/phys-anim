#! /usr/bin/env python

import os
import sys
import math
import string
import numpy as np


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
  s += '</box>\n';
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
  s += '</sphere>\n';
  return s;


if (__name__ == '__main__'):
  nx = 20;
  ny = 20;

  pos = np.empty((nx*ny, 3));
  xr = np.arange(nx) * 3.0 - nx*3.0/2.0;
  yr = np.arange(ny) * 4.0 + 5.0;
  ind = 0;
  for y in yr:
    for x in xr:
      pos[ind] = [x, y, 0];
      ind += 1;

  spheres = [];
  boxes = [];
  for p in pos:
    if (np.random.rand() < 0.5):
      spheres.append(p);
    else:
      boxes.append(p);

  materials = ['blue plastic', 'dark gray', 'light gray', 'astroglide', 'red', 'green', 'blue', 'dirt'];
  materials_xml = '  <materials>\
    <material name="blue plastic" density="1" friction="0.6" restitution="0.9" cr="0.5" cg="0.6" cb="0.9" />\
    <material name="dark gray" density="1" friction="0.9" restitution="0.9" cr="0.2" cg="0.2" cb="0.2" />\
    <material name="light gray" density="1" friction="0.9" restitution="0.9" cr="0.7" cg="0.7" cb="0.7" />\
    <material name="astroglide" density="1" friction="0.0" restitution="1.0" cr="0.2" cg="0.8" cb="0.4" />\
    <material name="red" density="1" friction="0.0" restitution="1.0" cr="0.8" cg="0.2" cb="0.4" />\
    <material name="green" density="1" friction="0.0" restitution="1.0" cr="0.2" cg="0.8" cb="0.4" />\
    <material name="blue" density="1" friction="0.0" restitution="1.0" cr="0.2" cg="0.4" cb="0.8" />\
    <material name="dirt"       density="0" friction="0.9" restitution="0.4" cr="0.8" cg="0.6" cb="0.4" />\
  </materials>\n';
  ground = '  <ground>\
    <bodymaterial name="dirt" />\
  </ground>\n'

  world = '<world>\n';
  world += materials_xml;
  world += '  <bodies>\n';
  world += ground;
  for s in spheres:
    world += sphere(s);
  for b in boxes:
    world += box(b);

  world += '  </bodies>\n';
  world += '</world>\n';

  print(world)

