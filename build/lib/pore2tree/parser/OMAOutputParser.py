import glob
import os
import tqdm

OUTPUT_FOLDER = 'Output'
OG_OUTPUT_FOLDER = os.join.path(OUTPUT_FOLDER, 'OrthologousGroupsFasta')

class OMAOutputError(Exception):
    pass

def add_arguments(arg_parser):
    '''
        Adds arguments to overall program's parser.
    '''
    arg_parser.add_argument('--standalone_path', required=True,
                            help='Path to the OMA standalone directory, which '
                                 'has completed.')

def process_arguments(args):
    '''
        Post-processes arguments to the program.
    '''
    setattr(args, 'oxml',
            os.path.join(args.standalone_path,
                         'Output',
                         'HierarchicalGroups.orthoxml'))

    # Add output_path
    setattr(args, 'output_path',
            os.path.join(args.standalone_path,
                         'Output',
                         'HOGPROP'))


class OMAOutputParser(object):

    def __init__(self, args):
        '''
        Initialise the OMAOutputParser
        '''
        self.args = args


    def __call__(self):
        '''
        This yields all the OGs, based on their respective output folder
        :return: 
        '''

        output_fn = os.path.join(self.args.standalone_path, OG_OUTPUT_FOLDER)

        for og in tqdm(glob.glob(output_fn+'*.fa'), desc='Loading OGs', unit=' og'):


