import importlib
from abc import ABC, abstractmethod
from textwrap import dedent
from collections.abc import Sequence

from astropy.table import QTable, Row, Table, TableAttribute

from sunpy.util.util import get_width

__all__ = ['BaseQueryResponse', 'QueryResponseTable', 'BaseClient']


class BaseQueryResponse(Sequence):
    """
    An Abstract Base Class for results returned from BaseClient.

    Notes
    -----
    * A QueryResponse object must be able to be instantiated with only one
      iterable argument. (i.e. the ``__init__`` must only have one required
      argument).
    * The `client` property must be settable.
    * The base class does not prescribe how you store the results from your
      client, only that it must be possible to represent them as an astropy
      table in the ``build_table`` method.
    * `__getitem__` **must** return an instance of the type it was called on.
      I.e. it must always return an object of ``type(self)``.

    """

    @abstractmethod
    def build_table(self):
        """
        Return an `astropy.table.Table` representation of the query response.
        """

    @property
    @abstractmethod
    def client(self):
        """
        An instance of `BaseClient` used to generate the results.

        Generally this is used to fetch the results later.

        .. note::

            In general, this doesn't have to be the same instance of
            ``BaseClient``, this is left to the client developer. If there is a
            significant connection overhead in creating an instance of a client
            you might want it to be the same instance as used for the search.
        """

    @client.setter
    @abstractmethod
    def client(self, value):
        pass

    @property
    @abstractmethod
    def blocks(self):
        """
        A `collections.abc.Sequence` object which contains the records
        contained within the Query Response.
        """

    def response_block_properties(self):
        """
        Returns a set of class attributes on all the response blocks.

        Returns
        -------
        s : `set`
            List of strings, containing attribute names in the response blocks.
        """
        set()

    def __str__(self):
        """Print out human-readable summary of records retrieved"""
        return '\n'.join(self.build_table().pformat(show_dtype=False))

    def __repr__(self):
        """Print out human-readable summary of records retrieved"""
        return object.__repr__(self) + "\n" + str(self)

    def _repr_html_(self):
        return self.build_table()._repr_html_()

    def show(self, *cols):
        """
        Returns response tables with desired columns for the Query.

        Parameters
        ----------
        \\*cols : `tuple`
            Name of columns to be shown.

        Returns
        -------
        `astropy.table.Table`
            A table showing values for specified columns.
        """
        table = self.build_table()
        if len(cols) == 0:
            return table
        tablecols = table.columns
        valid_cols = [col for col in cols if col in tablecols]
        return table[valid_cols]


class QueryResponseRow(Row):
    """
    A row subclass which knows about the client of the parent table.
    """

    def as_table(self):
        """
        Return this Row as a length one Table
        """
        return self.table[self.index:self.index + 1]


class QueryResponseTable(QTable, BaseQueryResponse):
    Row = QueryResponseRow

    client = TableAttribute()
    display_keys = TableAttribute(default=slice(None))
    hide_keys = TableAttribute()

    # This is a work around for https://github.com/astropy/astropy/pull/11217
    # Remove when min astropy version is > 4.2.1
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for attr in list(kwargs):
            descr = getattr(self.__class__, attr, None)
            if isinstance(descr, TableAttribute):
                setattr(self, attr, kwargs.pop(attr))

    def build_table(self):
        return self

    @property
    def blocks(self):
        return list(self.iterrows())

    def unhide_columns(self):
        self.display_keys = slice(None)
        self.hide_keys = None
        return self

    @property
    def _display_table(self):
        """
        Apply the display_keys and hide_keys attributes to the table.

        This removes any keys in hide keys and then slices by any keys in
        display_keys to return the correct table.
        """

        keys = list(self.colnames)
        if self.hide_keys:
            # Index only the keys not in hide keys in order
            [keys.remove(key) for key in self.hide_keys if key in keys]

        tslice = self.display_keys
        if tslice != slice(None):
            tslice = [dk for dk in self.display_keys if dk in keys]

        table = self[tslice]
        # The slicing operation resets display and hide keys to default, but we
        # have already applied it
        table.unhide_columns()

        return table

    def __str__(self):
        """Print out human-readable summary of records retrieved"""
        return '\n'.join(self._display_table.pformat(show_dtype=False))

    def __repr__(self):
        """Print out human-readable summary of records retrieved"""
        return object.__repr__(self) + "\n" + str(self._display_table)

    def _repr_html_(self):
        return QTable._repr_html_(self._display_table)

    def show(self, *cols):
        table = self.copy()
        table.unhide_columns()

        if len(cols) == 0:
            return table

        valid_cols = [col for col in cols if col in table.colnames]
        table = table[valid_cols]

        # The slicing operation resets display and hide keys to default, but we
        # want to bypass it here.
        table.unhide_columns()
        return table


def _print_client(client, html=False):
    """
    Given a BaseClient instance will print out each registered attribute.

    Parameters
    ----------
    client : `sunpy.net.base_client.BaseClient`
        The instance class to print for.
    html : bool
        Will return a html table instead.

    Returns
    -------
    `str`
        String with the client.
    """
    width = -1 if html else get_width()
    class_name = f"{client.__module__+'.' or ''}{client.__class__.__name__}"
    attrs = client.register_values()
    lines = []
    t = Table(names=["Attr Type", "Name", "Description"],
              dtype=["U80", "U80", "U80"])
    for client_key in attrs.keys():
        # Work around for * attrs having one length.
        if len(attrs[client_key]) == 1 and attrs[client_key][0] == "*":
            t.add_row((client_key.__name__, "All", "All valid values"))
            continue
        for name, desc in attrs[client_key]:
            t.add_row((client_key.__name__, name, desc))
    lines = [class_name, dedent(client.__doc__.partition("\n\n")[0])]
    if html:
        lines = [f"<p>{line}</p>" for line in lines]
    lines.extend(t.pformat_all(show_dtype=False, max_width=width, align="<", html=html))
    return '\n'.join(lines)


class BaseClient(ABC):
    """
    This defines the Abstract Base Class for each download client.

    The BaseClient has several abstract methods that ensure that any subclass enforces the bare minimum API.
    These are `search`, `fetch` and `_can_handle_query`.
    The last one ensures that each download client can be registered with Fido.

    Most download clients should subclass `~sunpy.net.dataretriever.GenericClient`.
    If the structure of `~sunpy.net.dataretriever.GenericClient`
    is not useful you should use `~sunpy.net.BaseClient`.
    `~sunpy.net.vso.VSOClient` and `~sunpy.net.jsoc.JSOCClient`
    are examples of download clients that subclass ``BaseClient``.
    """

    _registry = dict()

    def __init_subclass__(cls, *args, **kwargs):
        """
        An __init_subclass__ hook initializes all of the subclasses of a given class.
        So for each subclass, it will call this block of code on import.
        This replicates some metaclass magic without the need to be aware of metaclasses.
        Here we use this to register each subclass in a dict that has the `_can_handle_query` attribute.
        This is then passed into the UnifiedDownloaderFactory so we can register them.
        This means that Fido can use the clients internally.
        """
        super().__init_subclass__(**kwargs)

        # We do not want to register GenericClient since its a dummy client.
        if cls.__name__ in ('GenericClient'):
            return

        cls._registry[cls] = cls._can_handle_query

        if hasattr(cls, "_attrs_module"):
            from sunpy.net import attrs

            name, module = cls._attrs_module()
            module_obj = importlib.import_module(module)

            existing_mod = getattr(attrs, name, None)
            if existing_mod and existing_mod is not module_obj:
                raise NameError(f"{name} has already been registered as an attrs name.")

            setattr(attrs, name, module_obj)

            if name not in attrs.__all__:
                attrs.__all__.append(name)

        # Register client attrs after it has regsitered its own attrs
        from sunpy.net import attr
        values = cls.register_values()
        # If the client has no support, we won't try to register attrs
        if values:
            attr.Attr.update_values({cls: values})

    def __repr__(self):
        """
        Returns the normal repr plus the pretty client __str__.
        """
        return object.__repr__(self) + "\n" + str(self)

    def __str__(self):
        """
        This enables the "pretty" printing of BaseClient.
        """
        return _print_client(self)

    def _repr_html_(self):
        """
        This enables the "pretty" printing of the BaseClient with html.
        """
        return _print_client(self, html=True)

    @abstractmethod
    def search(self, *args, **kwargs):
        """
        This enables the user to search for data using the client.

        Must return a subclass of `BaseQueryResponse`.
        """

    @abstractmethod
    def fetch(self, *query_results, path=None, overwrite=False, progress=True,
              max_conn=5, downloader=None, wait=True, **kwargs):
        """
        This enables the user to fetch the data using the client, after a search.

        Parameters
        ----------
        query_results:
            Results to download.
        path : `str` or `pathlib.Path`, optional
            Path to the download directory
        overwrite : `bool`, optional
            Replace files with the same name if True.
        progress : `bool`, optional
            Print progress info to terminal.
        max_conns : `int`, optional
            Maximum number of download connections.
        downloader : `parfive.Downloader`, optional
            The download manager to use.
        wait : `bool`, optional
           If `False` ``downloader.download()`` will not be called. Only has
           any effect if `downloader` is not `None`.

        Returns
        -------
        `parfive.Results`
            The results object, can be `None` if ``wait`` is `False`.
        """

    @classmethod
    @abstractmethod
    def _can_handle_query(cls, *query):
        """
        This enables the client to register what kind of searches it can handle, to prevent Fido
        using the incorrect client.
        """

    @staticmethod
    def check_attr_types_in_query(query, required_attrs={}, optional_attrs={}):
        """
        Check a query againsted required and optional attributes.

        Returns `True` if *query* contains all the attrs in *required_attrs*,
        and if *query* contains only attrs in both *required_attrs* and *optional_attrs*.
        """
        query_attrs = {type(x) for x in query}
        all_attrs = required_attrs.union(optional_attrs)

        return required_attrs.issubset(query_attrs) and query_attrs.issubset(all_attrs)

    @classmethod
    def register_values(cls, *query):
        """
        This enables the client to register what kind of Attrs it can use directly.

        Returns
        -------
        `dict`
            A dictionary with key values of Attrs and the values are a tuple of
            ("Attr Type", "Name", "Description").
        """
        return {}
