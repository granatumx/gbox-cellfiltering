import matplotlib.pyplot as plt

from granatum_sdk import Granatum
import time

def main():
    tic = time.perf_counter()

    gn = Granatum()

    df = gn.pandas_from_assay(gn.get_import('assay'))
    mingenes = gn.get_arg('min_genes_per_cell')
    maxgenes = gn.get_arg('max_genes_per_cell')
    mt_percent = gn.get_arg('mt_genes_percent')

    uniquegenecount = df.astype(bool).sum(axis=0)
    totalgenecount = df.sum(axis=0)
    mtrows = df[df.index.str.startswith('MT')]
    mtgenecount = mtrows.sum(axis=0)
    mtpercent = mtgenecount.div(totalgenecount)
    colsmatching = uniquegenecount.T[(uniquegenecount.T >= mingenes) & (uniquegenecount.T <= maxgenes) & (mtpercent.T <= mt_percent)].index.values
    adata = df.loc[:, colsmatching]

    num_orig_cells = uniquegenecount.T.index.size
    num_filtered_cells = len(colsmatching)

    num_lt_min = uniquegenecount.T[(uniquegenecount.T < mingenes)].index.size
    num_gt_max = uniquegenecount.T[(uniquegenecount.T > maxgenes)].index.size
    num_gt_mt = uniquegenecount.T[(mtpercent.T > mt_percent)].index.size

    gn.add_result("Number of cells is now {} out of {} original cells with {} below min genes, {} above max genes, and {} above mt percentage threshold.".format(num_filtered_cells, num_orig_cells, num_lt_min, num_gt_max, num_gt_mt), "markdown")

    gn.export(gn.assay_from_pandas(adata), "Filtered Cells Assay", dynamic=False)

    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    timing = "* Finished cell filtering step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == '__main__':
    main()
