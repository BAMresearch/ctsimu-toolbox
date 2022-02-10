class Step:
    """ An image processing step to be run in the pipeline. """

    def __init__(self, stepIdentifier):
        self.identifier = None   # General description of operation
        self.prepared   = False  # Is step prepared to process first image?
        self.pipe       = None   # Pointer to the processing pipeline to which this step belongs.
        self.setIdentifier(stepIdentifier)

    def setIdentifier(self, identifier):
        self.identifier = identifier

    def setPipeline(self, pipePtr):
        self.pipe = pipePtr

    def getIdentifier(self):
        return self.identifier

    def setPrepared(self, prepared):
        self.prepared = prepared

    def isPrepared(self):
        return self.prepared

    # Virtual methods, implemented in children:
    def prepare(self):
        """ Prepare step before processing an image. """
        # Should be called once by run().
        pass

    def run(self, image):
        """ Run step's specified operation on given image. """
        pass

    def followUp(self):
        """ Follow-up procedures when all images have been processed. """
        pass