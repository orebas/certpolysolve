import subprocess
import sys
import tempfile
from pathlib import Path

def main():
    executable=sys.argv[1]
    with tempfile.TemporaryDirectory() as tmp:
        path=Path(tmp)/'input.pol'
        good=['Secular;\nRational;\nReal;\nDegree=2;\n2/3 0\n-2/3 3\n',
              'Secular; Degree=0;',
              'Secular; Real; Degree=1; 0.25 1.5']
        long_integer = '0' * 300000 + '1'
        literals = {
            'Integer': ['+1', '-2', long_integer],
            'Rational': ['+1/-2', '-3/+4', '1', '1/+' + long_integer],
            'FloatingPoint': ['.25', '1.', '-0.25', '+1e+2', '1.E-2', long_integer + '.0'],
        }
        for kind, values in literals.items():
            good.extend(f'Secular; Real; {kind}; Degree=1; {value} 0' for value in values)
        for text in good:
            path.write_text(text)
            r=subprocess.run([executable,'--secular',str(path)],capture_output=True,text=True,timeout=10)
            assert r.returncode==0 and 'goal_met=1' in r.stdout,(r.returncode,r.stdout,r.stderr)
        bad=['Degree=1; 1 0 0 0','Secular; Degree=2049;',
             'Secular; Real; Degree=1; 0 1',
             'Secular; Real; Degree=1; 1 0 extra',
             'Secular; Degree=1; Precision=15; 1 0 0 0',
             'Secular; Real; Degree=2; 1 0 2 0',
             'Secular; Integer; Real; Degree=1; 1.5 0',
             'Secular; Sparse; Degree=1; 1 0 0 0']
        invalid_literals = {
            'Integer': ['+', '--1', '1e0', '1.0', long_integer + 'x'],
            'Rational': ['1/', '/1', '1/0', '1/2/3', '1.0', '1/+' + long_integer + 'x'],
            'FloatingPoint': ['.', '-.', '1e', '1e+', '1..0', '1/2', '1@2', long_integer + '.0e-'],
        }
        for kind, values in invalid_literals.items():
            bad.extend(f'Secular; Real; {kind}; Degree=1; {value} 0' for value in values)
        for text in bad:
            path.write_text(text)
            r=subprocess.run([executable,'--secular',str(path)],capture_output=True,text=True,timeout=10)
            assert r.returncode==1,(text[:120],r.returncode,r.stdout,r.stderr)
    print('secular CLI checks passed')
if __name__=='__main__': main()
