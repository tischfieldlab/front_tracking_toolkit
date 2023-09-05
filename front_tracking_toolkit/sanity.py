import pandas as pd


def check_image_metadata_within_sample(img_meta: pd.DataFrame):
    ignore = ['AcquiredDate', 'CreationDate', 'index']
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
                    print(f'In sample {group}, key {key} changed from "{last[key]}" at index {last["index"]} to "{row[key]}" at index {row["index"]}')

            last = row.to_dict()


def check_image_metadata_across_samples(img_meta: pd.DataFrame):
    ignore = ['AcquiredDate', 'CreationDate', 'index', 'XMetresPerPixel', 'YMetresPerPixel']
    for stain, rows in img_meta.groupby('stain'):
        for c in rows.columns:
            if c not in ignore:
                uniq = rows[c].unique()
                print(c, len(uniq), uniq)

