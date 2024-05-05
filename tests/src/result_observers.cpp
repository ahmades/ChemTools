#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include <cassert>

#include <fmt/format.h>
#include <fmt/printf.h>

#include <catch2/catch.hpp>

#include <boost/filesystem.hpp>
#include <boost/system/error_code.hpp>

#include "results/console_logger.hpp"
#include "results/csv_writer.hpp"
#include "results/hdf5_writer.hpp"
#include "test_config.hpp"

// concrete subjects

class Motion : public Results::Subject
{
public:
  Motion(double const acceleration_, double const target_position_)
      : acceleration{acceleration_}, target_position{target_position_}, time{0.0}, position{0.0}, velocity{0.0}
  {
    RegisterResults();
  }

  void Compute()
  {
    time += 1.0;
    velocity = acceleration * time;
    position = velocity * time;
    distance_to_target = target_position - position;
    // notify
    Results::Subject::NotifyObserver();
  }

private:
  double acceleration;
  double target_position;
  double time;
  double position;
  double velocity;
  double distance_to_target;

  void RegisterResults()
  {
    // file name
    Results::Subject::results_meta.SetFileName("motion");

    // results
    {
      std::string const group("Main results");
      Results::Subject::results_meta.Register(group, "Time", "Time", "t", "s", &time);
      Results::Subject::results_meta.Register(group, "Position", "Position", "x", "m", &position);
      Results::Subject::results_meta.Register(group, "Velocity", "Velocity", "V", "m/s", &velocity);
    }

    // other results
    Results::Subject::results_meta.Register("Auxiliary results", "DistToTarg", "Distance to target", "Dx", "m", &distance_to_target);
  }
};

TEST_CASE("Results can be observed", "[results]")
{

  std::string const path = "./";

  Motion motion(0.25 // acceleration
                ,
                200.0); // target position
  size_t attached_observers = 0;
  {
    attached_observers = motion.AttachedObservers();
    REQUIRE(0 == attached_observers);
    if (TestConfig::verbose)
      fmt::print("Motion attached observer(s) = {} \n", attached_observers);

    std::unique_ptr<Results::CSVWriter> csv_writer;
    std::unique_ptr<Results::HDF5Writer> hdf5_writer;
    std::unique_ptr<Results::ConsoleLogger> console_logger;

    try
      {
        // attach to csv writer
        csv_writer = std::make_unique<Results::CSVWriter>(motion, path);
        attached_observers = motion.AttachedObservers();
        REQUIRE(1 == attached_observers);
        if (TestConfig::verbose)
          fmt::print("Motion attached observer(s) = {} \n", attached_observers);

        // attach to hdf5 writer
        hdf5_writer = std::make_unique<Results::HDF5Writer>(motion, path);
        hdf5_writer->SetCompressionLevel(4);
        attached_observers = motion.AttachedObservers();
        REQUIRE(2 == attached_observers);
        if (TestConfig::verbose)
          fmt::print("Motion attached observer(s) = {} \n", attached_observers);

        // attach to console logger
        console_logger = std::make_unique<Results::ConsoleLogger>(motion);
        attached_observers = motion.AttachedObservers();
        REQUIRE(3 == attached_observers);
        if (TestConfig::verbose)
          fmt::print("Motion attached observer(s) = {} \n", attached_observers);
      }
    catch (boost::filesystem::filesystem_error const& fs_err)
      {
        fmt::print("{}\n", fs_err.what());
        std::exit(EXIT_FAILURE);
      }

    size_t const n_iter = 30;     // number of iterations
    size_t const i_detach = 9;    // index of iteration where observer is detached
    size_t const i_reattach = 19; // index of iteration where observer is reattached

    for (size_t i = 0; i < n_iter; i++)
      {
        if (i == i_detach)
          {
            // detach from csv writer
            motion.DetachObserver(*csv_writer.get());
            attached_observers = motion.AttachedObservers();
            REQUIRE(2 == attached_observers);
            if (TestConfig::verbose)
              fmt::print("Motion attached observer(s) = {} \n", attached_observers);
          }
        else if (i == i_reattach)
          {
            // reattach from csv writer
            motion.AttachObserver(*csv_writer.get());
            attached_observers = motion.AttachedObservers();
            REQUIRE(3 == attached_observers);
            if (TestConfig::verbose)
              fmt::print("Motion attached observer(s) = {} \n", attached_observers);
          }
        motion.Compute();
      }
  }

  attached_observers = motion.AttachedObservers();
  REQUIRE(0 == attached_observers);
  if (TestConfig::verbose)
    fmt::print("Motion attached observer(s) = {} \n", attached_observers);
}
