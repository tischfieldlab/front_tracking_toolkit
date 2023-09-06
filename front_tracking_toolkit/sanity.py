import pandas as pd


def check_image_metadata_within_sample(img_meta: pd.DataFrame):
    # fields which we expect to change from one image to the next, within an image series
    ignore = [
        'AcquiredDate',
        'CreationDate',
        'index',
        'filename',
        'AcquiredDelta'
    ]
    print("I will check each sample for internal consistancy...")
    for group, rows in img_meta.groupby(['subject', 'stain']):
        last = None
        for _, row in rows.iterrows():
            if last is None:
                last = row.to_dict()
                continue

            for key in last.keys():
                if key in ignore:
                    continue

                if last[key] != row[key]:
                    print(f' -> WARN: In sample {group}, key "{key}" changed from "{last[key]}" at index {last["index"]} to "{row[key]}" at index {row["index"]}')

            last = row.to_dict()


def check_image_metadata_across_samples(img_meta: pd.DataFrame):
    # fields which we expect to change from one image to the next, within and across image series
    ignore = [
        'AcquiredDate',
        'CreationDate',
        'index',
        'XMetresPerPixel',
        'YMetresPerPixel',
        'filename',
        'AcquiredDelta',
        'subject',
    ]
    print("I will check across samples for consistancy...")
    for stain, rows in img_meta.groupby('stain'):
        print(f'Checking stain "{stain}..."')
        for c in rows.columns:
            if c not in ignore:
                uniq = rows[c].unique()
                if len(uniq) > 1:
                    print(f' -> WARN: Found {len(uniq)} distinct values for "{c}": {uniq}')

