class ProcessingStep:
    """ An image processing step to be run in the pipeline. """

    def __init__(self, stepIdentifier):
        self._identifier = None   # General description of operation
        self._prepared   = False  # Is step prepared to process first image?
        self._pipe       = None   # Pointer to the processing pipeline to which this step belongs.
        self.setIdentifier(stepIdentifier)

    def setIdentifier(self, identifier):
        self._identifier = identifier

    def setPipeline(self, pipePtr):
        self._pipe = pipePtr

    def getIdentifier(self):
        return self._identifier

    def setPrepared(self, prepared):
        self._prepared = prepared

    def isPrepared(self):
        return self._prepared

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