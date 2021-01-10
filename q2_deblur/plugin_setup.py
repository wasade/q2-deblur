# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import qiime2.plugin
import importlib
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Sequence
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (Sequences,
                                           SequencesWithQuality,
                                           PairedEndSequencesWithQuality,
                                           JoinedSequencesWithQuality)

import q2_deblur

from q2_deblur._type import DeblurStats
from q2_deblur._format import DeblurStatsFmt, DeblurStatsDirFmt


citations = qiime2.plugin.Citations.load('citations.bib', package='q2_deblur')

plugin = qiime2.plugin.Plugin(
    name='deblur',
    version=q2_deblur.__version__,
    website='https://github.com/biocore/deblur',
    package='q2_deblur',
    citations=[citations['amir2017deblur']],
    description=('This QIIME 2 plugin wraps the Deblur software for '
                 'performing sequence quality control.'),
    short_description='Plugin for sequence quality control with Deblur.'
)


plugin.register_formats(DeblurStatsFmt, DeblurStatsDirFmt)
plugin.register_semantic_types(DeblurStats)
plugin.register_semantic_type_to_format(DeblurStats,
                                        artifact_format=DeblurStatsDirFmt)


_parameter_descriptions = {
    'mean_error': ("The mean per nucleotide error, used for original "
                   "sequence estimate."),
    'indel_prob': ('Insertion/deletion (indel) probability (same for N '
                   'indels).'),
    'indel_max': "Maximum number of insertion/deletions.",
    'trim_length': "Sequence trim length, specify -1 to disable trimming.",
    'left_trim_len': "Sequence trimming from the 5' end. A value of 0 will "
                     "disable this trim.",
    'min_reads': ("Retain only features appearing at least min_reads "
                  "times across all samples in the resulting feature "
                  "table."),
    'min_size': ("In each sample, discard all features with an abundance "
                 "less than min_size."),
    'jobs_to_start': "Number of jobs to start (if to run in parallel).",
    'hashed_feature_ids': "If true, hash the feature IDs.",
    'sample_stats': "If true, gather stats per sample."
}


_output_descriptions = {
    'table': 'The resulting denoised feature table.',
    'representative_sequences': 'The resulting feature sequences.',
    'stats': 'Per-sample stats if requested.'
}


_parameters = {
    'mean_error': qiime2.plugin.Float,
    'indel_prob': qiime2.plugin.Float,
    'indel_max': qiime2.plugin.Int,
    'trim_length': qiime2.plugin.Int,
    'left_trim_len': qiime2.plugin.Int % qiime2.plugin.Range(0, None),
    'min_reads': qiime2.plugin.Int,
    'min_size': qiime2.plugin.Int,
    'jobs_to_start': qiime2.plugin.Int,
    'hashed_feature_ids': qiime2.plugin.Bool,
    'sample_stats': qiime2.plugin.Bool
}


_outputs = [('table', FeatureTable[Frequency]),
            ('representative_sequences', FeatureData[Sequence]),
            ('stats', DeblurStats)]


plugin.methods.register_function(
    function=q2_deblur.denoise_16S,
    inputs={
        'demultiplexed_seqs': SampleData[Sequences |
                                         SequencesWithQuality |
                                         PairedEndSequencesWithQuality |
                                         JoinedSequencesWithQuality],
    },
    parameters=_parameters,
    outputs=_outputs,
    input_descriptions={
        'demultiplexed_seqs': 'The demultiplexed sequences to be denoised.',
    },
    parameter_descriptions=_parameter_descriptions,
    output_descriptions=_output_descriptions,
    name='Deblur sequences using a 16S positive filter.',
    description=('Perform sequence quality control for Illumina data using '
                 'the Deblur workflow with a 16S reference as a positive '
                 'filter. Only forward reads are supported at this time. The '
                 'specific reference used is the 88% OTUs from Greengenes '
                 '13_8. This mode of operation should only be used when data '
                 'were generated from a 16S amplicon protocol on an Illumina '
                 'platform. The reference is only used to assess whether each '
                 'sequence is likely to be 16S by a local alignment using '
                 'SortMeRNA with a permissive e-value; the reference is not '
                 'used to characterize the sequences.')
)


plugin.methods.register_function(
    function=q2_deblur.denoise_other,
    inputs={
        'demultiplexed_seqs': SampleData[Sequences |
                                         SequencesWithQuality |
                                         PairedEndSequencesWithQuality |
                                         JoinedSequencesWithQuality],
        'reference_seqs': FeatureData[Sequence],
    },
    parameters=_parameters,
    outputs=_outputs,
    input_descriptions={
        'demultiplexed_seqs': 'The demultiplexed sequences to be denoised.',
        'reference_seqs': ("Positive filtering database. Keep all "
                           "sequences aligning to these sequences."),
    },
    parameter_descriptions=_parameter_descriptions,
    output_descriptions=_output_descriptions,
    name='Deblur sequences using a user-specified positive filter.',
    description=('Perform sequence quality control for Illumina data using '
                 'the Deblur workflow, including positive alignment-based '
                 'filtering. Only forward reads are supported at this time. '
                 'This mode of execution is particularly useful when '
                 'operating on non-16S data. For example, to apply Deblur to '
                 '18S data, you would want to specify a reference composed of '
                 '18S sequences in order to filter out sequences which do not '
                 'appear to be 18S. The assessment is performed by local '
                 'alignment using SortMeRNA with a permissive e-value '
                 'threshold.')
)

plugin.visualizers.register_function(
    function=q2_deblur.visualize_stats,
    inputs={'deblur_stats': DeblurStats},
    parameters={},
    input_descriptions={
        'deblur_stats': 'Summary statistics of the Deblur process.'
    },
    parameter_descriptions={},
    name='Visualize Deblur stats per sample.',
    description='Display Deblur statistics per sample'
)

importlib.import_module('q2_deblur._transformer')
