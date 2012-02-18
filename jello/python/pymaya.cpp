// TODO: change list arguments/returns to numpy nparrays

// boost python interface headers
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/extract.hpp>

#include <vector>
#include <string>

using namespace boost::python;


int check_length(list values) {
  int nvalue = len(values);
  return nvalue;
}

char const* hello_world() {
 return("hello, world");
}

BOOST_PYTHON_MODULE(pymaya) {
  using namespace boost::python;

  def("hello_world", hello_world);
  def("check_length", check_length);
}

