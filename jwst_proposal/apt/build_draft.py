"""Build a MIRI-only APT draft; no submission or target audit.

The schema representation follows a public MIRI/MRS APT file (program 1523).
No proposal-specific settings or cached timings are copied from that program.
APT must open the generated file and compute its own timing model.
"""
import argparse
import csv
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
NS = 'http://www.stsci.edu/JWST/APT'
MRS = NS + '/Template/MiriMRS'
XSI = 'http://www.w3.org/2001/XMLSchema-instance'
ET.register_namespace('', NS)
ET.register_namespace('mmrs', MRS)
ET.register_namespace('xsi', XSI)


def add(parent, name, text=None, ns=NS, **attrs):
    e = ET.SubElement(parent, '{' + ns + '}' + name, attrs)
    if text is not None:
        e.text = str(text)
    return e


def dec_string(degrees):
    sign = '+' if degrees >= 0 else '-'
    seconds = round(abs(degrees) * 3600, 3)
    d, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f'{sign}{int(d):02d} {int(m):02d} {s:06.3f}'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=HERE / 'work/generated',
                        help='Keep generated input separate from the APT-processed review file.')
    output = parser.parse_args().output_dir
    output.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader((BASE / 'inputs/jwst_sample_cycle6.csv').open()))
    sky_file = HERE / 'sky_fields/selected_sky_fields.csv'
    sky_fields = {r['id']: r for r in csv.DictReader(sky_file.open())} if sky_file.exists() else {}
    recipes = {r['id']: r for r in csv.DictReader((BASE / 'review/acquisition_recipes.csv').open())}
    assert len(recipes) == 24
    assert len(rows) == 24
    assert sum(r['id'].startswith('F') for r in rows) == 12
    assert sum(r['id'].startswith('R') for r in rows) == 12
    root = ET.Element('{' + NS + '}JwstProposal', schemaVersion='63',
                      APTVersion='Version 2024.1  JWST PRD: PRDOPSSOC-065, Roman PRD: RPRDDEVSOC-014',
                      PRDVersion='PRDOPSSOC-065')
    info = add(root, 'ProposalInformation')
    title_abstract = (BASE / 'apt_title_abstract.txt').read_text()
    add(info, 'Title', title_abstract.split('Title:\n', 1)[1].split('\n\nAbstract:', 1)[0].strip())
    abstract = title_abstract.split('Abstract:\n', 1)[1].strip()
    add(info, 'Abstract', abstract)
    # Do not write ProposalID. APT treats any value, including 0, as an already
    # assigned number, switches the Submission tool to resubmission mode and then
    # rejects it ('Proposal ID: 0 is too small'). A new proposal omits the element;
    # APT obtains the real number from STScI at first Submit (2026-09-29 fix).
    add(info, 'ProposalCategory', 'GO')
    add(info, 'ScientificCategory', 'Galaxies and the Intergalactic Medium')
    add(info, 'Cycle', '6')
    # APT's importer requires this container even for an anonymous timing draft.
    add(add(info, 'PrincipalInvestigator'), 'InvestigatorAddress')
    targets = add(root, 'Targets')
    requests = add(root, 'DataRequests')
    group = add(requests, 'ObservationGroup')
    add(group, 'Label', '24 AGN - MIRI/MRS')
    add(group, 'Comments', '12 fading and 12 rising AGNs. A/B/C groups 27/60/27, except R06 B=80; four dithers and one integration. Flux-dependent F560W/FND acquisitions follow review/acquisition_recipes.csv. Sky coordinates come from the reviewed sky-field table when available; otherwise 60-arcsec north offsets require review.')
    links = add(root, 'LinkingRequirements')
    inventory = []
    for idx, row in enumerate(rows):
        recipe = recipes[row['id']]
        assert recipe['target'] == row['target'], 'Acquisition recipe belongs to a different target'
        if row['id'] in sky_fields:
            assert sky_fields[row['id']]['target'] == row['target'], 'Sky field belongs to a different target'
        groups_a = 27
        groups_b = 80 if row['id'] == 'R06' else 60
        for sky in (False, True):
            num = 2 * idx + 1 + int(sky)
            label = row['id'] + ('-SKY' if sky else '')
            target = add(targets, 'Target', **{'{' + XSI + '}type': 'FixedTargetType'})
            add(target, 'Number', num)
            add(target, 'TargetName', label)
            if not sky:
                add(target, 'TargetArchiveName', row['target'])
            add(target, 'TargetID', label)
            sky_comment = 'Sky field selected from archival infrared catalogs and images; see sky_fields/sky_field_review.pdf.' if row['id'] in sky_fields else 'UNVETTED sky placeholder 60 arcsec north; choose clean field before submission.'
            add(target, 'Comments', sky_comment if sky else f"{row['target']}; {'fading' if row['id'].startswith('F') else 'rising'} AGN; source identity and coordinates from working sample.")
            add(target, 'ProperMotionNotApplicable', 'true')
            add(target, 'Extended', 'YES' if sky else 'NO')
            add(target, 'Category', 'Calibration' if sky else 'Galaxy')
            coords = row['ra_hms'].replace(':', ' ') + ' ' + (dec_string(float(row['dec']) + 1 / 60) if sky else row['dec_dms'].replace(':', ' '))
            if sky and row['id'] in sky_fields:
                coords = sky_fields[row['id']]['coordinates']
            coord_element = add(target, 'EquatorialCoordinates', Value=coords)
            for component in ['RAUncertainty', 'DecUncertainty']:
                uncertainty = add(coord_element, component)
                add(uncertainty, 'Value')
                add(uncertainty, 'Units', 'Arcsec')
            add(target, 'BackgroundTargetReq', 'false' if sky else 'true')
            if not sky:
                add(target, 'BackgroundTargets', f'{num+1} {label}-SKY')
            add(target, 'TargetConfirmationRun', 'false')
            obs = add(group, 'Observation', AutoTarget='false')
            add(obs, 'Number', num)
            add(obs, 'TargetID', f'{num} {label}')
            add(obs, 'Label', label + (' background' if sky else ' science'))
            add(obs, 'Instrument', 'MIRI')
            template = add(add(obs, 'Template'), 'MiriMRS', ns=MRS)
            def m(parent, key, value=None):
                return add(parent, key, value, ns=MRS)
            m(template, 'AcqTargetID', 'NONE' if sky else 'Same Target as Observation')
            if not sky:
                m(template, 'AcqFilter', recipe['filter'])
                m(template, 'AcqReadoutPattern', recipe['readout'])
                m(template, 'AcqGroups', recipe['groups'])
            m(template, 'Detector', 'MRS')
            dither = m(m(template, 'Dithers'), 'MrsDitherSpecification')
            m(dither, 'DitherType', '4-Point')
            m(dither, 'OptimizedFor', 'EXTENDED SOURCE' if sky else 'POINT SOURCE')
            m(dither, 'Direction', 'NEGATIVE')
            m(template, 'SimultaneousImaging', 'NO')
            m(template, 'Subarray', 'FULL')
            m(template, 'PrimaryChannel', 'ALL')
            m(template, 'GratingWheelADirection', 'NEUTRAL')
            exposures = m(template, 'ExposureList')
            for setting, groups in [('SHORT(A)', groups_a), ('MEDIUM(B)', groups_b), ('LONG(C)', 27)]:
                exposure = m(exposures, 'Exposure')
                for key, value in [('Exposures', 1), ('Wavelength', setting), ('Dither', 'Dither 1'),
                                   ('ReadoutPatternLong', 'FASTR1'), ('ReadoutPatternShort', 'FASTR1'),
                                   ('GroupsLong', groups), ('GroupsShort', groups),
                                   ('IntegrationsLong', 1), ('IntegrationsShort', 1)]:
                    m(exposure, key, value)
            add(obs, 'CoordinatedParallel', 'false')
            add(obs, 'SpecialRequirements')
            add(obs, 'Visit', Number='1')
            inventory.append(dict(observation=num, target=label, role='sky' if sky else 'science',
                                  coordinates=coords, extended='YES' if sky else 'NO',
                                  dither_optimized_for='EXTENDED SOURCE' if sky else 'POINT SOURCE',
                                  groups_A=groups_a, groups_B=groups_b, groups_C=27,
                                  integrations=1, dithers=4, target_acquisition='NONE' if sky else f"{recipe['filter']} {recipe['readout']} {recipe['groups']} groups"))
        link = add(links, 'GroupWithinLink', Number=str(idx + 1), Sequence='true', NonInterruptible='true')
        add(link, 'Observation', f"{row['id']} science (Obs {2*idx+1})")
        add(link, 'Observation', f"{row['id']}-SKY background (Obs {2*idx+2})")
        add(link, 'Time')
    ET.indent(root)
    xml_name = 'clagn24_miri_draft.xml'
    xml = ET.tostring(root, encoding='utf-8', xml_declaration=True)
    (output / xml_name).write_bytes(xml)
    with zipfile.ZipFile(output / 'clagn24_miri_draft.aptx', 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('manifest', f'\r\nName: AptAttributes\r\nProposalFilename: {xml_name}\r\n\r\n')
        z.writestr(xml_name, xml)
    with (output / 'observation_inventory.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(inventory[0]))
        writer.writeheader()
        writer.writerows(inventory)
    print('Wrote 24 AGN + 24 sky observations, linked in 24 pairs; APT timing not yet computed.')


if __name__ == '__main__':
    main()
