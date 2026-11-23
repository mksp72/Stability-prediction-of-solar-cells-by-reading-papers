"""Decode complete stability trajectories stored in one Excel cell.

Long JSON traces are compressed losslessly to fit Excel's 32,767-character
cell limit; ordinary shorter traces remain plain JSON.
"""
import base64,json,zlib

def decode_stability_data(value):
    if not value:return []
    if value.startswith('zlib64:'):
        value=zlib.decompress(base64.b64decode(value[len('zlib64:'):])).decode('utf8')
    return json.loads(value)
