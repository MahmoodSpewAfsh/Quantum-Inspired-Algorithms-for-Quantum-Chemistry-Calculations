class QAIA:
    r"""
    The base class of QAIA.

    This class contains the basic and common functions of all the algorithms.

    Args:
        J (Union[numpy.array, csr_matrix]): The coupling matrix with shape :math:`(N x N)`.
        h (numpy.array): The external field with shape :math:`(N x 1)`.
        x (numpy.array): The initialized spin value with shape :math:`(N x batch_size)`. Default: ``None``.
        n_iter (int): The number of iterations. Default: ``1000``.
        batch_size (int): The number of sampling. Default: ``1``.
    """

    # pylint: disable=too-many-arguments
    def __init__(self, J, h=None, x=None, n_iter=1000, batch_size=1):
        """Construct a QAIA algorithm."""
        self.J = J
        if h is not None and len(h.shape) < 2:
            h = h[:, np.newaxis]
        self.h = h
        self.x = x
        # The number of spins
        self.N = self.J.shape[0]
        self.n_iter = n_iter
        self.batch_size = batch_size

    def initialize(self):
        """Randomly initialize spin values."""
        self.x = 0.02 * (np.random.rand(self.N, self.batch_size) - 0.5)

    def calc_cut(self, x=None):
        r"""
        Calculate cut value.

        Args:
            x (numpy.array): The spin value with shape :math:`(N x batch_size)`.
                If ``None``, the initial spin will be used. Default: ``None``.
        """
        if x is None:
            sign = np.sign(self.x)
        else:
            sign = np.sign(x)

        return 0.25 * np.sum(self.J.dot(sign) * sign, axis=0) - 0.25 * self.J.sum()

    def calc_energy(self, x=None):
        r"""
        Calculate energy.

        Args:
            x (numpy.array): The spin value with shape :math:`(N x batch_size)`.
                If ``None``, the initial spin will be used. Default: ``None``.
        """
        if x is None:
            sign = np.sign(self.x)
        else:
            sign = np.sign(x)

        if self.h is None:
            return -0.5 * np.sum(self.J.dot(sign) * sign, axis=0), sign
        return -0.5 * np.sum(self.J.dot(sign) * sign, axis=0, keepdims=True) - self.h.T.dot(sign), sign


class OverflowException(Exception):
    r"""
    Custom exception class for handling overflow errors in numerical calculations.

    Args:
        message: Exception message string, defaults to "Overflow error".
    """

    def __init__(self, message="Overflow error"):
        self.message = message
        super().__init__(self.message)
"""Coherent Ising Machine with chaotic feedback control algorithm."""
# pylint: disable=invalid-name



class SFC(QAIA):
    r"""
    Coherent Ising Machine with separated feedback control algorithm.

    Reference: `Coherent Ising machines with optical error correction
    circuits <https://onlinelibrary.wiley.com/doi/full/10.1002/qute.202100077>`_.

    Args:
        J (Union[numpy.array, csr_matrix]): The coupling matrix with shape :math:`(N x N)`.
        h (numpy.array): The external field with shape :math:`(N, )`.
        x (numpy.array): The initialized spin value with shape :math:`(N x batch_size)`. Default: ``None``.
        n_iter (int): The number of iterations. Default: ``1000``.
        batch_size (int): The number of sampling. Default: ``1``.
        dt (float): The step size. Default: ``0.1``.
        k (float): parameter of deviation between mean-field and error variables. Default: ``0.2``.
    """

    # pylint: disable=too-many-arguments,too-many-instance-attributes
    def __init__(
        self,
        J,
        h=None,
        x=None,
        n_iter=1000,
        batch_size=1,
        dt=0.1,
        k=0.2,
    ):
        """Construct SFC algorithm."""
        super().__init__(J, h, x, n_iter, batch_size)
        self.J = csr_matrix(self.J)
        self.N = self.J.shape[0]
        self.dt = dt
        self.n_iter = n_iter
        self.k = k
        # pumping parameters
        self.p = np.linspace(-1, 1, self.n_iter)
        # coupling strength
        self.xi = np.sqrt(2 * self.N / np.sum(self.J**2))
        # rate of change of error variables
        self.beta = np.linspace(0.3, 0, self.n_iter)
        # coefficient of mean-field term
        self.c = np.linspace(1, 3, self.n_iter)
        self.initialize()

    def initialize(self):
        """Initialize spin values and error variables."""
        if self.x is None:
            self.x = np.random.normal(0, 0.1, (self.N, self.batch_size))
        self.e = np.zeros_like(self.x)

        if self.x.shape[0] != self.N:
            raise ValueError(f"The size of x {self.x.shape[0]} is not equal to the number of spins {self.N}")

    # pylint: disable=attribute-defined-outside-init
    def update(self):
        """Dynamical evolution."""
        for i in range(self.n_iter):
            if self.h is None:
                z = -self.xi * (self.J @ self.x)
            else:
                z = -self.xi * (self.J @ self.x + self.h)
            f = np.tanh(self.c[i] * z)
            self.x = self.x + (-self.x**3 + (self.p[i] - 1) * self.x - f - self.k * (z - self.e)) * self.dt
            self.e = self.e + (-self.beta[i] * (self.e - z)) * self.dt

            if np.isnan(self.x).any():
                raise OverflowException("Value is too large to handle due to large dt or xi.")
