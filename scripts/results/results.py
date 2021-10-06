import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt


def get_paths() :
    parser = argparse.ArgumentParser()
    parser.add_argument( "--file", "-f", type = str, required = True )
    parser.add_argument( "--work_dir", "-wd", type = str, required = True )
    args = parser.parse_args()
    return (args.file, args.work_dir)

class Result:
    def __init__( self, file_name, work_dir ) :
        self.work_dir = work_dir;
        self.h5_file = h5py.File( file_name, "r" )
    
    def print_keys(self) :
        print( "Keys: %s" % list( self.h5_file.keys() ) )
    
    def get_set( self, group_name, set_name ) :
        the_set = self.h5_file[group_name][set_name]
        attributes = {}
        for k in the_set.attrs.keys():
            attributes[k] = the_set.attrs[k]
        return attributes, np.array( the_set )

    def print_set( self, group_name, set_name ) :
        print( self.get_set( group_name, set_name ) )

    def plot( self,  x_group, x_set, y_group, y_sets, save_name ) :
        fig, ax = plt.subplots()
        x_attrs, x_data = self.get_set( x_group, x_set )
        for y_set in y_sets :
            y_attrs, y_data = self.get_set( y_group, y_set )
            ax.plot( x_data, y_data, label = y_attrs['Notation'] + ' (' + y_attrs["Unit"] + ')' )
        ax.set( xlabel = x_attrs["Notation"] + ' (' + x_attrs["Unit"] + ')' )
        if( len( y_sets ) == 1  ):
            ax.set( ylabel = y_attrs['Name'] + ' (' +  y_attrs['Unit'] + ')' )
        ax.legend()
        ax.grid()
        fig.savefig( self.work_dir + save_name )

