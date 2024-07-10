#ifndef TEST_CONFIG_HPP
#define TEST_CONFIG_HPP

namespace TestConfig
{

  // verbosity switch
#ifdef TEST_VERBOSE
  static bool const verbose = true;
#else
  static bool const verbose = false;
#endif

  // tolerance
  static double const tolerance = 100.0 * std::numeric_limits<float>::epsilon();

} // namespace TestConfig

#endif // TEST_CONFIG_HPP
