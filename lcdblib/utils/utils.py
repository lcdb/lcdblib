from collections.abc import Iterable

def flatten(iter):
    """
    Flatten an arbitrarily nested iterable whose innermost items are strings
    into a flat list of strings.
    """
    if isinstance(iter, dict):
        iter = iter.values()
    def gen():
        for item in iter:
            if isinstance(item, dict):
                item = item.values()
            if isinstance(item, Iterable) and not isinstance(item, str):
                yield from flatten(item)
            else:
                yield item
    return list(gen())


def test_flatten():
    assert sorted(flatten({
        'a': {
            'b': {
                'c': ['a','b','c'],
            },
        },
        'x': ['e', 'f', 'g'],
        'y': {
            'z': 'd'
        },
    })) == ['a', 'b', 'c', 'd', 'e', 'f', 'g']


def updatecopy(orig, update_with, keys=None, override=False):
    """
    Update a copy of a dictionary, with a bit more control than the built-in
    dict.update.

    Parameters
    -----------

    orig : dict
        Dict to update

    update_with : dict
        Dict with new values

    keys : list or None
        If not None, then only consider these keys in `update_with`. Otherwise
        consider all.

    override : bool
        If True, then this is similar to `dict.update`, except only those keys
        in `keys` will be considered. If False (default), then if a key exists
        in both `orig` and `update_with`, no updating will occur so `orig` will
        retain its original value.
    """
    d = orig.copy()
    if keys is None:
        keys = update_with.keys()
    for k in keys:
        if k in update_with:
            if k in d and not override:
                continue
            d[k] = update_with[k]
    return d
