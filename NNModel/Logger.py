# Code referenced from https://gist.github.com/gyglim/1f8dfb1b5c82627ae3efcfbbadb9f514
from torch.autograd import Variable
from torch.utils.tensorboard import SummaryWriter
import torch
import signal, os, sys
# import utils as u


def signal_term_handler(signal, frame):
    os.system("killall_tensorboard")
    print('got SIGTERM, killing tensorboard.')
    sys.exit(0)


try:
    from StringIO import StringIO  # Python 2.7
except ImportError:
    from io import BytesIO  # Python 3.x


def to_np(x):
    return x.data.cpu().numpy()


def to_var(x):
    if torch.cuda.is_available():
        x = x.cuda()
    return Variable(x)


class MetaLogger(object):
    """ This class encapsulates the tensorboard logging code, but uses the underlying Logger class
     to do all the real work."""

    def __init__(self, model, port=6006):
        self.model = model
        os.system("rm -rf " + os.getcwd() + "/logs")
        os.system("killall_tensorboard")
        os.system("tensorboard --logdir='" + os.getcwd() + "/logs' --port=" + str(port) + "&")
        self.writer = SummaryWriter(os.getcwd() + '/logs')

        # self.logger = Logger(os.getcwd()+'/logs')
        # kill previous opened process
        pid = os.system("lsof -i:6001 | grep tensorboa | awk '{print $2}'")
        if pid != None and pid > 0:
            os.system("kill -9 " + str(pid))

        os.system("firefox http://localhost:" + str(port) + "  &")

        signal.signal(signal.SIGTERM, signal_term_handler)

    def writeTensorboardLog(self, step, errTot, te, embeds):
        # step == epoch

        info = {"errTot": errTot}
        info["timePerEpoch"] = te

        for tag, value in info.items():
            # if tag.split(".")[0] != "BackHbond":
            #     continue
            try:
                for task, v in enumerate(value):
                    self.writer.add_scalar(tag + "_" + str(task), v, step + 1)
            except TypeError:
                self.writer.add_scalar(tag, value, step + 1)

        # visualize embeddings
        # self.writer.add_embedding(mat=embeds, metadata=list(u.resi_hash.values()), global_step=step)

        # (2) Log values and gradients of the parameters (histogram): here is the same as update_weights for "BackHbond"
        for tag, value in self.model.named_parameters():
            # print

            tag = tag.replace('.', '/')

            if type(value.grad) != type(None):
                self.writer.add_histogram(tag, to_np(value), step + 1)
                self.writer.add_histogram(tag + '/grad', to_np(value.grad), step + 1)
        return info

    def update_weights(self, step):
        for tag, value in self.model.named_parameters():
            # if tag.split(".")[0] != "BackHbond":
            #     continue

            tag = tag.replace('.', '/')

            if type(value.grad) != type(None):
                self.writer.add_histogram(tag, to_np(value), step + 1)
                self.writer.add_histogram(tag + '/grad', to_np(value.grad), step + 1)

    def shutdown(self):
        self.writer.flush()
        self.writer.close()
        print("Logger: killing tensorboard.")
        os.system("killall_tensorboard")


'''
class Logger(object):

	"""This class actually writes data on the tensorboard logger. It could be used
	directly or it can be used through the MetaLogger class"""

	def __init__(self, log_dir):
		"""Create a summary writer logging to log_dir."""
		#self.writer = tf.summary.FileWriter(log_dir)
		self.writer = tf.summary.create_file_writer(log_dir)
		#tf.logging.set_verbosity(tf.logging.WARN)
		os.environ['TF_CPP_MIN_LOG_LEVEL']='5'

	def list_summary(self, tags, values, step):
		assert len(tags) == len(values)
		for i, t in enumerate(tags):
			self.scalar_summary(t, values[i], step)

	def scalar_summary(self, tag, value, step):
		"""Log a scalar variable."""
		summary = tf.Summary(value=[tf.Summary.Value(tag=tag, simple_value=value)])
		self.writer.add_summary(summary, step)

	def image_summary(self, tag, images, step):
		"""Log a list of images."""

		img_summaries = []
		for i, img in enumerate(images):
			# Write the image to a string
			try:
				s = StringIO()
			except:
				s = BytesIO()
			scipy.misc.toimage(img).save(s, format="png")

			# Create an Image object
			img_sum = tf.Summary.Image(encoded_image_string=s.getvalue(),
									   height=img.shape[0],
									   width=img.shape[1])
			# Create a Summary value
			img_summaries.append(tf.Summary.Value(tag='%s/%d' % (tag, i), image=img_sum))

		# Create and write Summary
		summary = tf.Summary(value=img_summaries)
		self.writer.add_summary(summary, step)

	def histo_summary(self, tag, values, step, bins=1000):
		"""Log a histogram of the tensor of values."""

		# Create a histogram using numpy
		counts, bin_edges = np.histogram(values, bins=bins)

		# Fill the fields of the histogram proto
		hist = tf.HistogramProto()
		hist.min = float(np.min(values))
		hist.max = float(np.max(values))
		hist.num = int(np.prod(values.shape))
		hist.sum = float(np.sum(values))
		hist.sum_squares = float(np.sum(values**2))

		# Drop the start of the first bin
		bin_edges = bin_edges[1:]

		# Add bin edges and counts
		for edge in bin_edges:
			hist.bucket_limit.append(edge)
		for c in counts:
			hist.bucket.append(c)

		# Create and write Summary
		summary = tf.Summary(value=[tf.Summary.Value(tag=tag, histo=hist)])
		self.writer.add_summary(summary, step)
		self.writer.flush()
'''

