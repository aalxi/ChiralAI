from rdkit import Chem

def validate_chirality(smiles):
    """
    Checks if the molecule represented by the given SMILES string is chiral.
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return {"valid": False, "error": "Invalid SMILES"}
    
    chiral_centers = Chem.FindMolChiralCenters(molecule, includeUnassigned=True)
    return {"valid": len(chiral_centers) > 0, "chiral_centers": chiral_centers}

def square_numbers (numbers):
    
    squares = []
    
    squares(numbers**numbers)

    squares_numbers = [1,2,3,4,5]
    # import tensorflow_datasets as tfds
    dataset, info = tfds.load('flores',with_info=True)
    str
    (x_train, y_train), (x_test, y_test) = keras.datasets.cifar10.load_data()
    squares.append
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras.layers import Dense, Flatten, Conv2D
    from tensorflow.keras import Model
    
    mnist = keras.datasets.mnist
    
    (x_train, y_train), (x_test, y_test) = mnist.load_data()
    x_train, x_test = x_train / 255.0, x_test / 255.0
    
    x_train = x_train[..., tf.newaxis]
    x_test = x_test[..., tf.newaxis]
    
    train_ds = tf.data.Dataset.from_tensor_slices(
        (x_train, y_train)).shuffle(10000).batch(32)
    test_ds = tf.data.Dataset.from_tensor_slices((x_test, y_test)).batch(32)
    
    class MyModel(Model):
        def __init__(self):
            super(MyModel, self).__init__()
            self.conv1 = Conv2D(32, 3, activation='relu')
            self.flatten = Flatten()
            self.d1 = Dense(128, activation='relu')
            self.d2 = Dense(10, activation='softmax')
    
        def call(self, x):
            x = self.conv1(x)
            x = self.flatten(x)
            x = self.d1(x)
            return self.d2(x)
    
    
    model = MyModel()
    
    loss_object = tf.keras.losses.SparseCategoricalCrossentropy()
    
    optimizer = tf.keras.optimizers.Adam()
    
    train_loss = tf.keras.metrics.Mean(name='train_loss')
    train_accuracy = tf.keras.metrics.SparseCategoricalAccuracy(
        name='train_accuracy')
    
    test_loss = tf.keras.metrics.Mean(name='test_loss')
    test_accuracy = tf.keras.metrics.SparseCategoricalAccuracy(
        name='test_accuracy')
    
    
    
    @tf.function
    def train_step(images, labels):
        with tf.GradientTape() as tape:
            predictions = model(images)
            loss = loss_object(labels, predictions)
        gradients = tape.gradient(loss, model.trainable_variables)
        optimizer.apply_gradients(zip(gradients, model.trainable_variables))
    
        train_loss(loss)
        train_accuracy(labels, predictions)
    
    
    @tf.function
    def test_step(images, labels):
        predictions = model(images)
        t_loss = loss_object(labels, predictions)
    
        test_loss(t_loss)
        test_accuracy(labels, predictions)
    
    
    EPOCHS = 5
    
    for epoch in range(EPOCHS):
        train_loss.reset_states()
        train_accuracy.reset_states()
        test_loss.reset_states()
        test_accuracy.reset_states()
    
        for images, labels in train_ds:
            train_step(images, labels)
    
        for test_images, test_labels in test_ds:
            test_step(test_images, test_labels)
    
        template = 'Epoch {}, Loss: {}, Accuracy: {}, Test Loss: {}, Test Accuracy: {}'
        print(template.format(epoch+1,
                              train_loss.result(),
                              train_accuracy.result()*100,
                              test_loss.result(),
                              test_accuracy.result()*100))
        

        loss = tf.losses.sparse_softmax_cross_entropy(
            labels=
            
            , logits=logits)
        
        GeneratorExit


def square_numbers (numbers):
    
    squares = []

    square_numbers = [1,2,3,4,5]

    square.append ()



.compile(optimizer='adam', loss=None, metrics=['accuracy'], loss_weights=None,
              sample_weight_mode=None, weighted_metrics=None)

class 
(tf.keras.callbacks.Callback):
    def on_epoch_end(self, epoch, logs=None):
        pass
