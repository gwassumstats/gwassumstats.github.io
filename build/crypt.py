"""The encrypted part of the site's data: AES-256-GCM under a key derived from the passphrase with PBKDF2-HMAC-SHA256.

site/app.js derives the same key and decrypts with WebCrypto. The passphrase is short, so the many PBKDF2 rounds are
what slows down guessing it; the salt is new for every encryption, so a key is never used twice.
"""
import base64
import hashlib
import os

from cryptography.hazmat.primitives.ciphers.aead import AESGCM

ITERATIONS = 1_000_000  # about half a second per guess in a browser


def derive(passphrase: str, salt: bytes, iterations: int) -> bytes:
    return hashlib.pbkdf2_hmac("sha256", passphrase.encode(), salt, iterations, 32)


def encrypt(plain: bytes, passphrase: str, iterations: int = ITERATIONS) -> dict:
    """The blob as written to the site: parameters, salt, nonce, and ciphertext with its tag, in base64."""
    salt, iv = os.urandom(16), os.urandom(12)
    data = AESGCM(derive(passphrase, salt, iterations)).encrypt(iv, plain, None)
    return {"kdf": "PBKDF2-HMAC-SHA256", "iterations": iterations, "cipher": "AES-256-GCM",
            **{k: base64.b64encode(v).decode() for k, v in (("salt", salt), ("iv", iv), ("data", data))}}


def decrypt(blob: dict, passphrase: str) -> bytes:
    """Raises cryptography.exceptions.InvalidTag for a wrong passphrase or a changed blob."""
    salt, iv, data = (base64.b64decode(blob[k]) for k in ("salt", "iv", "data"))
    return AESGCM(derive(passphrase, salt, blob["iterations"])).decrypt(iv, data, None)
