#!/Users/homapapoli/miniconda3/envs/sortseq_analysis/bin/python3.10

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import mavenn

# Choose dataset
data_name = 'sortseq'

print(f"Loading dataset '{data_name}' ")

# Load dataset
data_df = load_example_dataset(data_name)

# Get and report sequence length
L = len(data_df.loc[0, 'x'])
print(f'Sequence length: {L:d} nucleotides')

# Split dataset
trainval_df, test_df = mavenn.split_dataset(data_df)

class AdditiveGPMapLayer(GPMapLayer):
    """Represents an additive G-P map."""

    @handle_errors
    def __init__(self, *args, **kwargs):
        """Construct layer instance."""

        # Call superclass constructor
        super().__init__(*args, **kwargs)

        """Build layer."""
        # Define theta_0
        self.theta_0 = self.add_weight(name='theta_0',
                                       shape=(1,),
                                       initializer=Constant(0.),
                                       trainable=True,
                                       regularizer=self.regularizer)

        # Define theta_lc parameters
        theta_lc_shape = (1, self.L, self.C)
        theta_lc_init = np.random.randn(*theta_lc_shape)/np.sqrt(self.L)
        self.theta_lc = self.add_weight(name='theta_lc',
                                        shape=theta_lc_shape,
                                        initializer=Constant(theta_lc_init),
                                        trainable=True,
                                        regularizer=self.regularizer)

    def call(self, x_lc):
        """Process layer input and return output."""
        # Shape input
        x_lc = tf.reshape(x_lc, [-1, self.L, self.C])

        phi = self.theta_0 + \
              tf.reshape(K.sum(self.theta_lc * x_lc, axis=[1, 2]),
                         shape=[-1, 1])

        return phi

