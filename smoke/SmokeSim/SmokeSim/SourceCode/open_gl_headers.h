// Created by Peter Kutz.
// You might need to change these includes anyway, depending on your setup.

#pragma once

#ifdef __APPLE__

// MAC:
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#else

// WINDOWS:
#include <windows.h>
#include <gl/GL.h>
#include <gl/GLU.h>
#include <gl/glut.h>

#endif
