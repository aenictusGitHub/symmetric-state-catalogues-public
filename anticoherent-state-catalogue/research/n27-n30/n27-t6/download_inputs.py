"""Download missing frozen certificate inputs using only the Python standard library."""
from pathlib import Path
import hashlib
import sys
import urllib.request

BASE_URL = "https://raw.githubusercontent.com/aenictusGitHub/symmetric-state-catalogues-public/main/anticoherent-state-catalogue/research/n27-n30/n27-t6/"
INPUT_SHA256 = {'MANIFEST.txt': 'abb5197bb7cabdb2f8a89d3d176134470272e04aca79c42ec0ba18d81d4f9656', 'README.txt': 'bd6c58a6b7f9dc108c0a097ff077d1e2c5c94ef674b1b2ceb5453f9e60e2a034', 'angular.py': '26bef4efb9711f8691775749b1babfdd0838ccafb116214c60d88263e5ec03fa', 'audit_ppt_recoupling.py': 'db57c56830ddbe33a5c24ea8e1a148c78db85d947840189779a9b620380bbec0', 'certificate.json': 'ca3b31ae137c1034edf1553a110fcdfb09b0b021d71a23333ac5acf2118e6f3e', 'cfp2.npz': 'd518a11a59829af211f4c2824de2d9d26be9a568eeb2357c133c75beb3c2e6f7', 'cfp3.npz': 'dc33bb6eab1896f5b924a6785502fb3d58099232563106515f38108748441d22', 'cfp4.npz': '25db221537e28ab45b102460b631bd070a92e29ccbfe1bdf88418068cab7206b', 'cfp5.npz': '91c05cdf3b7704f1c14b432fd38ea0908bc1ebd7c41e0a4c81052970faff3a10', 'cfp6-part1.npz': '112b9db010d168281fbe5df90920c5ef2626369e3a5b5610998374c7d4c8c6b7', 'cfp6-part2.npz': '18fe43038b48860d07ea8e875951e8fa87cfa01835baa2573670c45be1ac25f1', 'cfp6-part3.npz': 'bf029ef5606f7a886b529b60a62786390ee88c6c895bb673dcb4a1a9087496d2', 'cfp6-part4.npz': '7bc1d6d945e0b8fe96935f59cff96cc92941cfe7a844e15785d60a11a88fa3fe', 'cfp6-part5.npz': '240c41ab552bc2c8dc5aff8599f5728fca7ce3ee88eb6c437e9b0c55a6743b37', 'cfp6-part6.npz': 'c6108416475d344a6e974b0da6c8ea0f29e644aeb33db1de1c05f95e5db4f51f', 'congruences.npz': 'a8759b603d1019e1e5e52fc4a40a4935e30803dbd5f01be33319f057bec3a87b', 'dual-part1.npz': 'af2c6e43d548d321fe4283176101c7ee6c0cde786a800384033fbadb3ba0ef91', 'dual-part2.npz': '295f68ec61570d8fec5a0ca11cfb91a81672cbf7365dde1746b64e57fd5b0a31', 'ppt-factors.npz': 'c7584d94219f4f0956aac7f297b371942abcb7f5ff765512e1a417cfc7c23821', 'ppt_first_term.py': 'ec5a51d2ba4c568592cb08131384703b7111de2f08356abad3131be9dd18b74d', 'proof.tex': 'e92d39d1b0635d15a70d4ce0882ad9cc5933d610ae19c183140f88b9ba274bb2', 'recoupling-control.json': 'd17381a4737263892b9de4a4e3d2f46e91f9ebf5e22b578f25c30da8b8cc6dec', 'validated_bootstrap.py': '85fdcfe041b5df375288d875c3a18269740c4f301898e449e0cce823fdf8736f', 'verify_certificate.py': 'db3505bf2a5d6a578040f2dd6abfa49d2ca8d81efda126178c1067a829dfcf53'}

def main():
    root = Path(__file__).resolve().parent
    base = sys.argv[1].rstrip('/')+'/' if len(sys.argv)>1 else BASE_URL
    for name, expected in INPUT_SHA256.items():
        target = root/name
        if target.is_file() and hashlib.sha256(target.read_bytes()).hexdigest()==expected:
            print('Verified existing', name, flush=True)
            continue
        temporary = root/(name+'.download')
        digest = hashlib.sha256()
        try:
            request = urllib.request.Request(base+name, headers={'User-Agent':'AC-certificate-downloader/1.0'})
            with urllib.request.urlopen(request, timeout=180) as response, temporary.open('wb') as output:
                while True:
                    chunk = response.read(1024*1024)
                    if not chunk: break
                    output.write(chunk)
                    digest.update(chunk)
            if digest.hexdigest()!=expected:
                raise ValueError('SHA256 mismatch for '+name)
            temporary.replace(target)
            print('Downloaded and verified', name, flush=True)
        finally:
            if temporary.exists(): temporary.unlink()
    print('All immutable inputs verified. Run: python3 verify_certificate.py')

if __name__=='__main__':
    main()
