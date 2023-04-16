#!/usr/bin/env python

import sys
import os
import shutil
import unittest
import numpy as np
import matplotlib
import anndata
matplotlib.use("agg")

if os.environ.get("TEST_MODE", "INSTALL") == "DEV":
    sys.path.insert(0, "..")
import Cell_BLAST as cb
cb.config.RANDOM_SEED = 0

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"


class DirectiTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = anndata.read_h5ad("pollen.h5ad")
        cb.data.normalize(cls.data)

    def tearDown(self):
        if os.path.exists("./test_directi"):
            shutil.rmtree("./test_directi")

    def test_gau(self):
        model = cb.directi.fit_DIRECTi(
            self.data, latent_dim=10, epoch=3, path="./test_directi"
        )
        self.data.obsm["X_latent"] = model.inference(self.data)
        self.assertFalse(np.any(np.isnan(self.data.obsm["X_latent"])))
        with self.assertRaises(Exception):
            model.clustering(self.data)
        model.save()
        del model
        model = cb.directi.DIRECTi.load("./test_directi")
        latent2 = model.inference(self.data)
        self.assertTrue(np.all(self.data.obsm["X_latent"] == latent2))

    def test_catgau(self):
        model = cb.directi.fit_DIRECTi(
            self.data, latent_dim=10, cat_dim=20, epoch=3, path="./test_directi"
        )
        self.data.obsm["X_latent"] = model.inference(self.data)
        self.assertFalse(np.any(np.isnan(self.data.obsm["X_latent"])))
        _ = model.clustering(self.data)
        model.save()
        del model
        model = cb.directi.DIRECTi.load("./test_directi")
        latent2 = model.inference(self.data)
        self.assertTrue(np.all(self.data.obsm["X_latent"] == latent2))
        random_state = np.random.RandomState(0)
        model.gene_grad(
            self.data, random_state.randn(self.data.shape[0], 10)
        )

    def test_semisupervised_catgau(self):
        self.data.obs["cell_type1"] = self.data.obs["cell_type1"].astype(str)
        _ = cb.data.annotation_confidence(
            self.data, "cell_type1", used_vars=self.data.uns["scmap_genes"],
            return_group_percentile=False
        )
        self.data.obs.loc[
            cb.data.annotation_confidence(
                self.data, "cell_type1", return_group_percentile=True
            )[1] <= 0.5, "cell_type1"
        ] = ""
        model = cb.directi.fit_DIRECTi(
            self.data, genes=self.data.uns["scmap_genes"],
            latent_dim=10, cat_dim=20, prob_module="ZINB",
            supervision="cell_type1", epoch=3, path="./test_directi"
        )
        self.data.obsm["X_latent"] = model.inference(self.data)
        self.assertFalse(np.any(np.isnan(self.data.obsm["X_latent"])))
        model.save()
        del model
        model = cb.directi.DIRECTi.load("./test_directi")
        latent2 = model.inference(self.data)
        self.assertTrue(np.all(self.data.obsm["X_latent"] == latent2))

    def test_rmbatch(self):
        model = cb.directi.fit_DIRECTi(
            self.data, genes=self.data.uns["scmap_genes"],
            latent_dim=10, batch_effect="cell_type1",  # Just for test
            epoch=3, path="./test_directi"
        )
        self.data.obsm["X_latent"] = model.inference(self.data)
        self.assertFalse(np.any(np.isnan(self.data.obsm["X_latent"])))
        model.save()
        del model
        model = cb.directi.DIRECTi.load("./test_directi")
        latent2 = model.inference(self.data)
        self.assertTrue(np.all(self.data.obsm["X_latent"] == latent2))

        model = cb.directi.fit_DIRECTi(
            self.data, genes=self.data.uns["scmap_genes"],
            prob_module="NB", latent_dim=10, depth=0,
            batch_effect="cell_type1",  # Just for test
            rmbatch_module="RMBatch",
            epoch=3, path="./test_directi"
        )
        self.data.obsm["X_latent"] = model.inference(self.data)
        self.assertFalse(np.any(np.isnan(self.data.obsm["X_latent"])))
        model.save()
        del model
        model = cb.directi.DIRECTi.load("./test_directi")
        latent2 = model.inference(self.data)
        self.assertTrue(np.all(self.data.obsm["X_latent"] == latent2))

        model = cb.directi.fit_DIRECTi(
            self.data, genes=self.data.uns["scmap_genes"],
            prob_module="ZINB",
            latent_dim=10, batch_effect="cell_type1",  # Just for test
            rmbatch_module="MNN",
            epoch=3, path="./test_directi"
        )
        self.data.obsm["X_latent"] = model.inference(self.data)
        self.assertFalse(np.any(np.isnan(self.data.obsm["X_latent"])))
        model.save()
        del model
        model = cb.directi.DIRECTi.load("./test_directi")
        latent2 = model.inference(self.data)
        self.assertTrue(np.all(self.data.obsm["X_latent"] == latent2))

        model = cb.directi.fit_DIRECTi(
            self.data, genes=self.data.uns["scmap_genes"],
            prob_module="ZILN",
            latent_dim=10, batch_effect="cell_type1",  # Just for test
            rmbatch_module="MNNAdversarial",
            epoch=3, path="./test_directi"
        )
        self.data.obsm["X_latent"] = model.inference(self.data)
        self.assertFalse(np.any(np.isnan(self.data.obsm["X_latent"])))
        model.save()
        del model
        model = cb.directi.DIRECTi.load("./test_directi")
        latent2 = model.inference(self.data)
        self.assertTrue(np.all(self.data.obsm["X_latent"] == latent2))

        model = cb.directi.fit_DIRECTi(
            self.data, genes=self.data.uns["scmap_genes"],
            prob_module="LN",
            latent_dim=10, batch_effect="cell_type1",  # Just for test
            rmbatch_module="AdaptiveMNNAdversarial",
            epoch=3, path="./test_directi"
        )
        self.data.obsm["X_latent"] = model.inference(self.data)
        self.assertFalse(np.any(np.isnan(self.data.obsm["X_latent"])))
        model.save()
        del model
        model = cb.directi.DIRECTi.load("./test_directi")
        latent2 = model.inference(self.data)
        self.assertTrue(np.all(self.data.obsm["X_latent"] == latent2))


if __name__ == "__main__":
    # DirectiTest.setUpClass()
    # test = DirectiTest()
    # test.test_semisupervised_catgau()
    # test.tearDown()
    unittest.main()
