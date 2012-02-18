%module example
%{
extern int cube(int n);
%}
%include example.c
