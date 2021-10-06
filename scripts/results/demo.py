from results import *
       
if __name__ == '__main__':
    file_name, work_dir = get_paths()
    result = Result( file_name, work_dir )
    result.print_keys()
    result.print_set( 'Temperature', 'Temperature' )
    result.plot( 'Time', 'Time', 'Mass fractions', ('H2', 'O2', 'H2O'), 'mass_fractions.png' )
    result.plot( 'Time', 'Time', 'Temperature', ('Temperature',), 'temperature.png' )    
