from vac2fost.vac2fost import MCFOSTUtils


def test_unicity():
    found = []
    for block, lines in MCFOSTUtils.blocks_descriptors.items():
        for line in lines:
            print('from test : ', block, line)
            found += [param for param in line]

    for p in found:
        if found.count(p) > 1:
            print(f'non unique key detected {p}')
    assert len(set(found)) == len(found)

def test_writter_null():
    MCFOSTUtils.write_mcfost_conf('tests/output/writter_out.para')

def test_writter_args():
    MCFOSTUtils.write_mcfost_conf(
        'tests/output/writter_out_2.para',
        custom={'nphot_sed': 2}
    )

