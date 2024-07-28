#include <array>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include <fmt/format.h>
#include <fmt/printf.h>

#include <boost/filesystem.hpp>
#include <boost/system/error_code.hpp>

#include "results/hdf5_writer.hpp"

namespace Results
{

  HDF5Writer::HDF5Writer(Subject& subject_,
                         std::string const& path_)
      : subject(subject_),
        path(path_),
        file(nullptr),
        map_scalars(),
        compression_level{3}
  {
    // attach
    subject.AttachObserver(*this);

    try
      {
        Initialise();
      }
    catch (boost::filesystem::filesystem_error const& fs_err)
      {
        fmt::print("{}\n", fs_err.what());
        std::exit(EXIT_FAILURE);
      }
    catch (H5::FileIException const& error)
      { // catch failure caused by the H5File operations
        error.printErrorStack();
        std::exit(EXIT_FAILURE);
      }
    catch (H5::DataSetIException const& error)
      { // catch failure caused by the DataSet operations
        error.printErrorStack();
        std::exit(EXIT_FAILURE);
      }
    catch (H5::DataSpaceIException const& error)
      { // catch failure caused by the DataSpace operations
        error.printErrorStack();
        std::exit(EXIT_FAILURE);
      }
    catch (H5::AttributeIException const& error)
      { // catch failure caused by the Attribute operations
        error.printErrorStack();
        std::exit(EXIT_FAILURE);
      }
  }

  HDF5Writer::~HDF5Writer()
  {
    subject.DetachObserver(*this);
  }

  // Modifies the default compression level
  void HDF5Writer::SetCompressionLevel(int const compression_level_)
  {
    compression_level = compression_level_;
  }

  static void AppendScalar(H5::DataSet* const dataset,
                           double* const values)
  {
    // dataspace
    int const rank = 1;
    hsize_t dims[rank] = {1};
    hsize_t max_dims[rank] = {H5S_UNLIMITED};
    H5::DataSpace const mem_space(rank, dims, max_dims);

    // get the current number of elements in the dataset
    hsize_t const current_n_elems = dataset->getSpace().getSimpleExtentNpoints();

    // extend the dataset
    hsize_t const new_n_elems[rank] = {current_n_elems + 1};
    dataset->extend(new_n_elems);

    // select hyperslab
    H5::DataSpace const& file_space = dataset->getSpace();
    hsize_t const offset[rank] = {current_n_elems};
    file_space.selectHyperslab(H5S_SELECT_SET,
                               dims,
                               offset);

    // write
    dataset->write(&*values,
                   H5::PredType::NATIVE_DOUBLE,
                   mem_space,
                   file_space);
  }

  // updates the values of registered variables
  void HDF5Writer::Update(Subject& a_subject)
  {
    if (&a_subject == &subject)
      {
        for (auto it = map_scalars.begin(); it != map_scalars.end(); ++it)
          {
            AppendScalar(it->first.get(), it->second);
          }
      }
  }

  // bundles the different initialisation steps
  void HDF5Writer::Initialise()
  {
    OpenFile();
    CreateGroups();
    CreateScalarDataSets();
  }

  // opens the results file
  void HDF5Writer::OpenFile()
  {
    // check path
    boost::filesystem::path const fs_path(path);
    boost::system::error_code error_code;
    if (!boost::filesystem::exists(fs_path, error_code))
      {
        BOOST_FILESYSTEM_THROW(boost::filesystem::filesystem_error("HDF5Writer error",
                                                                   fs_path,
                                                                   error_code));
      }

    // still need to get absolute path if relative is supplied

    // create the file given the path and file name
    file = std::make_unique<H5::H5File>(path + '/' + subject.results_meta.file_name + ".h5", H5F_ACC_TRUNC);
  }

  void HDF5Writer::CreateGroups()
  {
    std::vector<Meta::Scalar> const& scalars = subject.results_meta.scalars;
    // extract unique group names and and create groups
    std::unordered_set<std::string> unique_groups;
    // scalar
    for (auto it_scalars = scalars.begin(); it_scalars != scalars.end(); ++it_scalars)
      {
        std::string const& group_name = std::get<Meta::group>(*it_scalars);
        if (unique_groups.find(group_name) == unique_groups.end())
          {
            unique_groups.insert(group_name);
            file->createGroup('/' + group_name);
          }
      }
  }

  void HDF5Writer::CreateScalarDataSets()
  {
    int const rank = 1;
    hsize_t dims[rank] = {0};
    hsize_t chunk_dims[rank] = {1000};
    hsize_t max_dims[rank] = {H5S_UNLIMITED};

    std::vector<Meta::Scalar> const& scalars = subject.results_meta.scalars;

    for (auto it_scalars = scalars.begin(); it_scalars != scalars.end(); ++it_scalars)
      {
        // retrieve meta
        std::string const& group = std::get<Meta::group>(*it_scalars);
        std::string const& set = std::get<Meta::set>(*it_scalars);
        std::string const& name = std::get<Meta::name>(*it_scalars);
        std::string const& notation = std::get<Meta::notation>(*it_scalars);
        std::string const& unit = std::get<Meta::unit>(*it_scalars);
        double* const value = std::get<Meta::value>(*it_scalars);

        // create data space
        H5::DataSpace dataspace(rank, dims, max_dims);

        // create property list
        H5::DSetCreatPropList prop_list;
        prop_list.setChunk(rank, chunk_dims);
        prop_list.setDeflate(6); // compression level, default is 0 = uncompressed.

        // create dataset
        std::unique_ptr<H5::DataSet> dataset =
            std::make_unique<H5::DataSet>(file->createDataSet(group + '/' + set,
                                                              H5::PredType::NATIVE_DOUBLE,
                                                              dataspace,
                                                              prop_list));

        // add attributes to dataset
        {
          // construct map of attributes
          std::unordered_map<std::string, std::string> map = {
              {"Name", name},
              {"Notation", notation},
              {"Unit", unit}};

          // create variable length string datatype
          H5::StrType const str_type(0, H5T_VARIABLE);
          H5::DataSpace const dataspace(H5S_SCALAR);

          for (auto& element : map)
            {
              // create arttribute using map key
              H5::Attribute attribute = dataset->createAttribute(element.first,
                                                                 str_type,
                                                                 dataspace);

              // write attribute value using map value
              attribute.write(str_type, element.second);
            }
        }

        // insert dataset and value in the scalars map
        map_scalars.insert(std::make_pair(std::move(dataset), value));
      }
  }

} // namespace Results
