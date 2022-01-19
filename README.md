# CoSAP
Comparative Sequencing Analysis Platform

This repository hosts the development of the CoSAP library.

## About CoSAP

CoSAP is an easy yet comprehensive pipeline creation tool for NGS data. It aims to provide a deeper insight about the powers and limitations of the current tools by
allowing users to compare results of different pipelines. 


## Requirements and Installation

CoSAP requires [Miniconda](https://docs.conda.io/en/latest/miniconda.html) to install required packages.
Onec you install the miniconda,
clone the CoSAP repository:
```bash
git@github.com:MBaysanLab/cosap.git
```
cd into CoSAP:
```bash
cd cosap
```
and run: 
```bash
make install
```

---

## First contact with Keras

The core data structures of Keras are __layers__ and __models__.
The simplest type of model is the [`Sequential` model](/guides/sequential_model/), a linear stack of layers.
For more complex architectures, you should use the [Keras functional API](/guides/functional_api/),
which allows to build arbitrary graphs of layers, or [write models entirely from scratch via subclasssing](/guides/making_new_layers_and_models_via_subclassing/).

Here is the `Sequential` model:

```python
from tensorflow.keras.models import Sequential

model = Sequential()
```

Stacking layers is as easy as `.add()`:

```python
from tensorflow.keras.layers import Dense

model.add(Dense(units=64, activation='relu'))
model.add(Dense(units=10, activation='softmax'))
```

Once your model looks good, configure its learning process with `.compile()`:

```python
model.compile(loss='categorical_crossentropy',
              optimizer='sgd',
              metrics=['accuracy'])
```

If you need to, you can further configure your optimizer. The Keras philosophy is to keep simple things simple,
while allowing the user to be fully in control when they need to (the ultimate control being the easy extensibility of the source code via subclassing).

```python
model.compile(loss=tf.keras.losses.categorical_crossentropy,
              optimizer=tf.keras.optimizers.SGD(
                  learning_rate=0.01, momentum=0.9, nesterov=True))
```

You can now iterate on your training data in batches:

```python
# x_train and y_train are Numpy arrays.
model.fit(x_train, y_train, epochs=5, batch_size=32)
```

Evaluate your test loss and metrics in one line:

```python
loss_and_metrics = model.evaluate(x_test, y_test, batch_size=128)
```

Or generate predictions on new data:

```python
classes = model.predict(x_test, batch_size=128)
```

What you just saw is the most elementary way to use Keras.

However, Keras is also a highly-flexible framework suitable to iterate on state-of-the-art research ideas.
Keras follows the principle of **progressive disclosure of complexity**: it makes it easy to get started,
yet it makes it possible to handle arbitrarily advanced use cases,
only requiring incremental learning at each step.

In much the same way that you were able to train & evaluate a simple neural network above in a few lines,
you can use Keras to quickly develop new training procedures or exotic model architectures.
Here's a low-level training loop example, combining Keras functionality with the TensorFlow `GradientTape`:

```python
import tensorflow as tf

# Prepare an optimizer.
optimizer = tf.keras.optimizers.Adam()
# Prepare a loss function.
loss_fn = tf.keras.losses.kl_divergence

# Iterate over the batches of a dataset.
for inputs, targets in dataset:
    # Open a GradientTape.
    with tf.GradientTape() as tape:
        # Forward pass.
        predictions = model(inputs)
        # Compute the loss value for this batch.
        loss_value = loss_fn(targets, predictions)

    # Get gradients of loss wrt the weights.
    gradients = tape.gradient(loss_value, model.trainable_weights)
    # Update the weights of the model.
    optimizer.apply_gradients(zip(gradients, model.trainable_weights))
```

For more in-depth tutorials about Keras, you can check out:

- [Introduction to Keras for engineers](https://keras.io/getting_started/intro_to_keras_for_engineers/)
- [Introduction to Keras for researchers](https://keras.io/getting_started/intro_to_keras_for_researchers/)
- [Developer guides](https://keras.io/guides/)

---

## Installation

Keras comes packaged with TensorFlow 2 as `tensorflow.keras`.
To start using Keras, simply [install TensorFlow 2](https://www.tensorflow.org/install).

---

## Release and compatibility

Keras has **nightly releases** (`keras-nightly` on PyPI)
and **stable releases** (`keras` on PyPI).
The nightly Keras releases are usually compatible with the corresponding version
of the `tf-nightly` releases
(e.g. `keras-nightly==2.7.0.dev2021100607` should be
used with `tf-nightly==2.7.0.dev2021100607`).
We don't maintain backward compatibility for nightly releases.
For stable releases, each Keras
version maps to a specific stable version of TensorFlow.

The table below shows the compatibility version mapping
between TensorFlow versions and Keras versions.

All the release branches can be found on [Github](https://github.com/keras-team/keras/releases).

All the release binaries can be found on [Pypi](https://pypi.org/project/keras/#history).

| Keras release | Note      | Compatible Tensorflow version |
| -----------   | ----------- | -----------        |
| [2.4](https://github.com/keras-team/keras/releases/tag/2.4.0)  | Last stable release of multi-backend Keras | < 2.5
| 2.5-pre| Pre-release (not formal) for standalone Keras repo | >= 2.5 < 2.6
| [2.6](https://github.com/keras-team/keras/releases/tag/v2.6.0)    | First formal release of standalone Keras.  | >= 2.6 < 2.7
| [2.7](https://github.com/keras-team/keras/releases/tag/v2.7.0-rc0)    | (Upcoming release) | >= 2.7 < 2.8
| nightly|                                            | tf-nightly

---
## Support

You can ask questions and join the development discussion:

- In the [TensorFlow forum](https://discuss.tensorflow.org/).
- On the [Keras Google group](https://groups.google.com/forum/#!forum/keras-users).
- On the [Keras Slack channel](https://kerasteam.slack.com). Use [this link](https://keras-slack-autojoin.herokuapp.com/) to request an invitation to the channel.

---

## Opening an issue

You can also post **bug reports and feature requests** (only)
in [GitHub issues](https://github.com/keras-team/keras/issues).


---

## Opening a PR

We welcome contributions! Before opening a PR, please read
[our contributor guide](https://github.com/keras-team/keras/blob/master/CONTRIBUTING.md),
and the [API design guideline](https://github.com/keras-team/governance/blob/master/keras_api_design_guidelines.md).
